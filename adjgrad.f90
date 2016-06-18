      SUBROUTINE ADJGRAD(MASTER,MYID,MYNID,MYCOMM, NSPACE,NDOF,       &   
                         LVERB,LCOVD,IFREQ,ISRC, AZMOD,FREQ,          &
                         PY, &
                         EPSS,DTMAX,FMAT_DIST, RCV,MSH,INV,MID,       &
                         GRADL,IERR)
!
!     Updates the gradient calculating: 
!       grad = grad - Re{ adj(J) delta d} = grad - Re{ adj(F) S^{-1} R^T delta d}
!     where R^T maps the observations to an [ndof x 1] vector.  The calculation 
!     follows three simple steps:
!
!       (1) Backpropagating the residuals: lambda = adj(S^{-1}) R^T (d - Ru) 
!       (2) Multiplies adj(F) lambda
!       (3) Take  Re{adj(F) S^{-1} R^T delta d} and stack
!
!     Realize the gradient points in the direction of max increase of the error 
!     function.  So it must be negated when calculation a search direction. 
!
!     The notation follows from Metivier et. al. 2012 but if one is more comfortable
!     with Pratt (1998) then think of lambda as the virtual force v.  The reason for 
!     the change in notation is I'd like to do a full Newton inversion at some point 
!     in time - B. Baker December 2012
! 
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'fwd_struc.h'
      TYPE (CMUMPS_STRUC) MID
      TYPE (RECV_INFO) RCV
      TYPE (INV_INFO) INV
      TYPE (MESH_INFO) MSH
      COMPLEX*8, INTENT(IN) :: FMAT_DIST(NSPACE)
      REAL*8, INTENT(IN) :: FREQ, AZMOD, PY
      REAL*4, INTENT(IN) :: EPSS, DTMAX
      INTEGER*4, INTENT(IN) :: MASTER,MYID,MYNID,MYCOMM,  &
                               NSPACE,NDOF, IFREQ,ISRC
      LOGICAL*4, INTENT(IN) :: LVERB, LCOVD
      INTEGER*4, INTENT(OUT) :: IERR
      REAL*4, INTENT(INOUT) :: GRADL(*)
!.... local variables
      COMPLEX*8, ALLOCATABLE :: ADJLAM(:)
      REAL*4, ALLOCATABLE :: GPROD(:) 
      LOGICAL*4 LDEBUG
      PARAMETER(LDEBUG = .FALSE.) 
!
!----------------------------------------------------------------------------------------!
!
!.... backpropgate
      IF (MYID == MASTER .AND. LVERB) WRITE(*,*) 'adjgrad: Backpropagating residuals...'
      ALLOCATE(ADJLAM(NDOF))
      CALL CADJLAM(MASTER,MYNID,MYCOMM, NDOF, LCOVD, IFREQ,ISRC, DTMAX,  &
                   AZMOD,EPSS,FREQ,PY, RCV,INV,MID, ADJLAM,IERR)
!
!.... multiply adj(F) v
      IF (MYID == MASTER .AND. LVERB) &
      WRITE(*,*) 'adjgrad: Correlating residuals with model derivatives...'
      IF (MYNID == MASTER) THEN
         ALLOCATE(GPROD(inv%NA35))
      ELSE
         ALLOCATE(GPROD(1))
      ENDIF
      CALL ADJFLAM_DIST(MASTER,MYCOMM, inv%NNPGL,NDOF,inv%NZ_FDIST,inv%NA35, inv%MYGRAD, &
                        inv%JCSC_FDIST,inv%ICSC_FDIST,FMAT_DIST, ADJLAM, GPROD) 

      IF (LDEBUG .AND. MYID == MASTER) THEN
         CALL PLOT_ADJGRAD_VTK(NGNOD,NDIM,msh%NEN, NDIM,msh%NNPG,msh%NDOF,       &
                               msh%NLXI,msh%NLETA, msh%NELEM,NDIM, ISRC,1,       &
                               FREQ, msh%LM,msh%IENG, msh%XLOCS,msh%ZLOCS, ADJLAM) 

      ENDIF
      DEALLOCATE(ADJLAM)
!
!.... stack to complete -[dS/dm u] S^{-1} R^T delta d
      IF (MYID == MASTER .AND. LVERB) WRITE(*,*) 'adjgrad: Updating gradient...'
      IF (MYNID == MASTER) CALL SAXPY(inv%NA35,-1.0,GPROD,1,GRADL,1) 
!     IF (MYNID == MASTER) THEN
!        DO 1 IA35=1,inv%NA35
!           GRADL(IA35) = GRADL(IA35) + REAL(CPROD(IA35))  
!   1    CONTINUE
!     ENDIF
      IF (ALLOCATED(GPROD)) DEALLOCATE(GPROD) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE ADJFLAM_DIST(MASTER,MYCOMM, NNPGL,NDOF,NZ_FDIST,NA35, MYGRAD,      &
                              JCSC_FDIST,ICSC_FDIST,FMAT_DIST, ADJLAM, GPROD)
!
!     Calculates the distributed matrix vector
!       adj(F) lambda = -adj[dS/dm_1 u, dS/dm_2 u, ..., dS/dm_m u] lambda 
!     where dS/dm is held in a distributed form
!
!     INPUT       MEANING
!     -----       ------- 
!     ADJLAM      solution of adjoint problem adj(S) lambda = R^T(d - Ru)
!     ICSC_FDIST  distributed CSC row pointer for FMAT_DIST 
!     JCSC_FDIST  distributed CSC column pointer for FMAT_DIST 
!     FMAT_DIST   local columns of F
!     MASTER      master process ID 
!     MYCOMM      MPI communicator 
!     MYGRAD      maps local columns to points in gradient
!     NA35        number of points in gradient
!     NNPGL       number of local points in gradient, local columns of F
!     NZ_FDIST    number of non-zeros in distributed F matrix
!
!     OUTPUT      MEANING
!     ------      ------- 
!     GPROD       product on master process
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      COMPLEX*8, INTENT(IN) :: FMAT_DIST(NZ_FDIST)
      COMPLEX*8 ADJLAM(NDOF) 
      INTEGER*4, INTENT(IN) :: MYGRAD(NNPGL), JCSC_FDIST(NNPGL+1), ICSC_FDIST(NZ_FDIST), &
                               NNPGL, NA35, NDOF, NZ_FDIST, MYCOMM, MASTER
      REAL*4, INTENT(OUT) :: GPROD(*)
!.... local variables
      REAL*4, ALLOCATABLE :: VWORK(:) 
      COMPLEX*8 CZERO
      INTEGER*4 MPIERR, I1,I2,I,INPINV,INPGL, JDOF
      PARAMETER(CZERO = CMPLX(0.0,0.0))
!
!----------------------------------------------------------------------------------------!
!
!.... sparse Hermitian matrix vector multiply
      ALLOCATE(VWORK(NA35))
      VWORK(1:NA35) = 0.0 !CZERO
      DO 1 INPGL=1,NNPGL
         INPINV = MYGRAD(INPGL)
         I1 = JCSC_FDIST(INPGL) 
         I2 = JCSC_FDIST(INPGL+1) - 1
         DO 2 I=I1,I2
            JDOF = ICSC_FDIST(I) 
            VWORK(INPINV) = VWORK(INPINV) - REAL(CONJG(FMAT_DIST(I))*ADJLAM(JDOF))
    2    CONTINUE 
    1 CONTINUE
!.... reduce onto master
      !CALL MPI_REDUCE(VWORK,CPROD,NA35,MPI_COMPLEX, MPI_SUM,MASTER, MYCOMM,MPIERR)  
      CALL MPI_REDUCE(VWORK,GPROD,NA35,MPI_REAL, MPI_SUM,MASTER, MYCOMM,MPIERR)
      DEALLOCATE(VWORK) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE CADJLAM(MASTER,MYNID,MYCOMM, NDOF, LCOVD, IFREQ,ISRC, DTMAX,  & 
                         AZMOD,EPSS,FREQ,PY, RCV,INV,MID, ADJLAM,IERR) 
!
!     This backpropagates the residuals.  This is the first step in ADJGRAD but to be 
!     clear when programming we instead just solve Metivier's 2.24
!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth (degrees)
!     FREQ       current frequency (Hz)
!     IFREQ      current frequency 
!     IRESTP     residual type (1) phase, (2) amplitude, (3) phase and amplitude
!     ISRC       source type
!     NDOF       number of equations in global matrix
!     MASTER     master process ID
!     MID        MUMPS structure
!     MYCOMM     MPI communicator
!     MYNID      process ID on communicator
!
!     OUTPUT     MEANING
!     ------     ------- 
!     ADJLAM     solution of adj(S) lambda = R^T(d - u)
!     IERR       error flag 
!
!.... variable declarations
      implicit none
      INCLUDE 'mpif.h'
!     INCLUDE 'mesh_inv.inc'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'fwd_struc.h'
      TYPE (CMUMPS_STRUC) MID
      TYPE (RECV_INFO) RCV
      TYPE (INV_INFO)  INV
      REAL*8, INTENT(IN) :: FREQ, AZMOD, PY
      REAL*4, INTENT(IN) :: EPSS, DTMAX
      LOGICAL*4, INTENT(IN) :: LCOVD
      INTEGER*4, INTENT(IN) :: MASTER,MYNID,MYCOMM, NDOF, IFREQ, ISRC 
      COMPLEX*8, INTENT(OUT) :: ADJLAM(NDOF) 
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      COMPLEX*8, ALLOCATABLE :: RESID(:,:), ESTR(:,:), OBS(:,:), RHS_DENSE(:) 
      REAL*8, ALLOCATABLE :: COVD8(:,:)
      INTEGER*4, ALLOCATABLE :: IPERM_REC(:)
      COMPLEX*8 CZERO, COSAZ, SINAZ, QN, QE, QZ, Q1, Q2, CPHM2CM, CFACT
      REAL*8 PI180
      REAL*4 UMAG, VMAG, WMAG, NMAG, EMAG, ZMAG, UPH, VPH, WPH, NPH, EPH, ZPH 
      INTEGER*4 NWORK,NRHSS,LDC,ISAVE9,ISAVE20,I,IREC,IWARN,MPIERR
      PARAMETER(CZERO = CMPLX(0.0,0.0)) 
      PARAMETER(PI180 = 0.017453292519943295D0)
      REAL*8, PARAMETER :: TWOPI = 6.2831853071795862
      LOGICAL*4, PARAMETER :: LSHIFT = .TRUE.
!
!----------------------------------------------------------------------------------------!
!
!.... innitialization and memory handling 
      IERR = 0
      IF (MYNID == MASTER) THEN
         NRHSS   = MID%NRHS      !number of RHSs
         ISAVE9  = MID%ICNTL(9)  !save transpose problem 
         ISAVE20 = MID%ICNTL(20) !sparse RHS?
         IF (ISAVE20 == 0) THEN 
            NWORK = NRHSS*MID%LRHS
            ALLOCATE(RHS_DENSE(NWORK))
            RHS_DENSE(1:NWORK) = MID%RHS(1:NWORK)
            DEALLOCATE(MID%RHS)
            ALLOCATE(MID%RHS(MID%N))
         ENDIF
         MID%ICNTL(9)  = 0 !want to solve transpose problem
         MID%ICNTL(20) = 0 !dense RHS
         MID%NRHS = 1      !only one RHS
         IF (LCOVD) THEN
            ALLOCATE(IPERM_REC(NDIM*rcv%NREC))
            LDC = 2*rcv%NREC*NDIM
            ALLOCATE(COVD8(LDC,LDC))
            COVD8(:,:) = 0.D0
            !CALL LAPLACE_1D(LDC,rcv%NREC,inv%DX, COVD8,IERR)
            CALL LAPLACE_1DV2(LDC,NDIM,rcv%NREC,inv%LDWGHT, inv%DX, &
                              inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC), COVD8,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'cadjlam: Error generating difference matrix!'
               RETURN
            ENDIF 
            CALL PERM_OBS(NDIM,NDIM,rcv%NREC, inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC), &
                          IPERM_REC,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'cadjlam: Error calling perm_obs!'
               RETURN
            ENDIF
         ENDIF
      ENDIF
!
!.... fill the RHS, this generates R^T delta d* 
      IF (MYNID == MASTER) THEN 
         ALLOCATE(RESID(NDIM,rcv%NREC)) 
         ALLOCATE(ESTR(NDIM,rcv%NREC))
         ALLOCATE(OBS(NDIM,rcv%NREC))
!
!....... put estimates in (N,E,Z) frame 
         COSAZ = CMPLX(SNGL(DCOS(AZMOD*PI180)),0.0)
         SINAZ = CMPLX(SNGL(DSIN(AZMOD*PI180)),0.0)
         DO 1 IREC=1,rcv%NREC
            IF (inv%LUNWRAP) THEN !need data, estimates as a complex number, est->(N,E,Z) 
               UMAG = REAL(inv%EST(1,IFREQ,IREC,ISRC))
               VMAG = REAL(inv%EST(2,IFREQ,IREC,ISRC)) 
               WMAG = REAL(inv%EST(3,IFREQ,IREC,ISRC)) 
               UPH  = IMAG(inv%EST(1,IFREQ,IREC,ISRC))
               VPH  = IMAG(inv%EST(2,IFREQ,IREC,ISRC))
               WPH  = IMAG(inv%EST(3,IFREQ,IREC,ISRC)) 
               !CALL CROTATE(THETA,Q1,Q2, RQ1,RQ2)
               Q1 = CPHM2CM(UMAG,UPH)
               Q2 = CPHM2CM(VMAG,VPH)
               CALL CROTATE(SNGL(AZMOD),Q1,Q2, ESTR(1,IREC),ESTR(2,IREC))
               !ESTR(1,IREC) = CPHM2CM(UMAG,UPH)*COSAZ - CPHM2CM(VMAG,VPH)*SINAZ
               !ESTR(2,IREC) = CPHM2CM(UMAG,UPH)*SINAZ + CPHM2CM(VMAG,VPH)*COSAZ  
               ESTR(3,IREC) = CPHM2CM(WMAG,WPH)  
               NMAG = REAL(inv%OBS(1,IFREQ,IREC,ISRC))
               EMAG = REAL(inv%OBS(2,IFREQ,IREC,ISRC))
               ZMAG = REAL(inv%OBS(3,IFREQ,IREC,ISRC))
               NPH  = IMAG(inv%OBS(1,IFREQ,IREC,ISRC))
               EPH  = IMAG(inv%OBS(2,IFREQ,IREC,ISRC))
               ZPH  = IMAG(inv%OBS(3,IFREQ,IREC,ISRC))
               OBS(1,IREC) = CPHM2CM(NMAG,NPH) 
               OBS(2,IREC) = CPHM2CM(EMAG,EPH)
               OBS(3,IREC) = CPHM2CM(ZMAG,ZPH)
            ELSE !data/est already complex just rotate est -> (N,E,Z)
               Q1 = inv%EST(1,IFREQ,IREC,ISRC)
               Q2 = inv%EST(2,IFREQ,IREC,ISRC)
               CALL CROTATE(SNGL(AZMOD),Q1,Q2, ESTR(1,IREC),ESTR(2,IREC)) 
               !ESTR(1,IREC) = inv%EST(1,IFREQ,IREC,ISRC)*COSAZ &
               !             - inv%EST(2,IFREQ,IREC,ISRC)*SINAZ 
               !ESTR(2,IREC) = inv%EST(1,IFREQ,IREC,ISRC)*SINAZ &
               !             + inv%EST(2,IFREQ,IREC,ISRC)*COSAZ 
               ESTR(3,IREC) = inv%EST(3,IFREQ,IREC,ISRC)  
               OBS(1,IREC) = inv%OBS(1,IFREQ,IREC,ISRC) 
               OBS(2,IREC) = inv%OBS(2,IFREQ,IREC,ISRC)
               OBS(3,IREC) = inv%OBS(3,IFREQ,IREC,ISRC)
            ENDIF
!
!.......... put data nd observations back into (u,v,w) frame completely
            IF (LSHIFT) THEN
               CFACT = CEXP(CMPLX(0.D0,-SNGL(TWOPI*FREQ*PY*rcv%YREC(IREC))))
                OBS(1:3,IREC) =  OBS(1:3,IREC)*CFACT
               ESTR(1:3,IREC) = ESTR(1:3,IREC)*CFACT
            ENDIF
    1    CONTINUE
!
!....... set residuals to backpropagate in (N,E,Z) frame
         IF (LCOVD) THEN
            CALL BPRESID_COVD4(NDIM,LDC, NDIM,rcv%NREC,inv%IRESTP,   IPERM_REC, &
                               DTMAX,FREQ,COVD8,inv%WGHTS(1:NDIM,IFREQ,1:rcv%NREC,ISRC), &
                               OBS,ESTR, RESID)
         ELSE
            CALL BPRESID4(NDIM,rcv%NREC,NDIM, inv%NORM,inv%IRESTP, EPSS, &  
                          FREQ,DTMAX, inv%WGHTS(1:NDIM,IFREQ,1:rcv%NREC,ISRC), &
                          OBS,ESTR, RESID)
         ENDIF
!
!....... unrotate residuals (u,v,w) frame so we are in same coordinates as S 
         DO 2 IREC=1,rcv%NREC
            QN = RESID(1,IREC)
            QE = RESID(2,IREC)
            QZ = RESID(3,IREC)
            CALL CROTATE(-SNGL(AZMOD),QN,QE, RESID(1,IREC),RESID(2,IREC))
            !RESID(1,IREC) = QN*COSAZ + QE*SINAZ
            !RESID(2,IREC) =-QN*SINAZ + QE*COSAZ
            RESID(3,IREC) = QZ 
    2    CONTINUE
!
!....... put into MUMPS RHS
         CALL RESID_DRHS(NDIM,NDOF,rcv%NREC,NDIM, inv%IBPHASE,FREQ,rcv%MRDOF, RESID,   &
                         mid%RHS,IWARN)
         IF (IWARN > 0)  &
         WRITE(*,*) 'cadjlam: You should check your DOF numbers for the receivers!'
         DEALLOCATE(RESID) 
         DEALLOCATE(ESTR) 
         DEALLOCATE(OBS) 
      ENDIF
!
!.... solution phase for S^T lambda* = R^T(d - R u)*  
      MID%JOB = 3 
      CALL CMUMPS(MID)
      IF (MID%INFO(1) < 0) THEN
         IF (MYNID == MASTER) WRITE(*,*) 'cadjlam: Error in solution phase'
         IERR = 1
         RETURN
      ENDIF
!
!.... have lambda* want lambda
      IF (MYNID == MASTER) THEN
         DO 3 I=1,MID%N
            ADJLAM(I) = CONJG(MID%RHS(I))
    3    CONTINUE
      ENDIF
      CALL MPI_BCAST(ADJLAM,NDOF,MPI_COMPLEX, MASTER,MYCOMM,MPIERR)
! 
!.... restore memory 
      IF (MYNID == MASTER) THEN
         MID%NRHS = NRHSS
         MID%ICNTL(9)  = ISAVE9
         MID%ICNTL(20) = ISAVE20
         DEALLOCATE(MID%RHS)
         IF (MID%ICNTL(20) == 0) THEN
            ALLOCATE(MID%RHS(NWORK))
            MID%RHS(1:NWORK) = RHS_DENSE(1:NWORK)
            DEALLOCATE(RHS_DENSE)
         ENDIF
         IF (ALLOCATED(IPERM_REC)) DEALLOCATE(IPERM_REC)
         IF (ALLOCATED(COVD8)) DEALLOCATE(COVD8)
         !IF (MYID == MASTER) WRITE(*,*)
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE RESID_RHS(MDIM,NDOF,NREC,NDIM, IBPHASE,FREQ,MRDOF, RESID,  &
                           IRHS_PTR,IRHS_SPARSE,RHS_SPARSE, IWARN)
! 
!     This routine calculate the residual populating the RHS for backpropgation.  Here, 
!     we use the principal of superposition to force our RHS into one vector.   
! 
!        x_1 = inv(A) b_1
!        x_2 = inv(A) b_2
!          . 
!          . 
!          .  
!        x_n = inv(A) b_n
! 
!     Then sum to obtain 
!        x_1 + x_2 + ... + x_n = inv(A) b_1 + inv(A) b_2 + ... inv(A) b_n 
!                              = inv(A) (b_1 + b_2 + ... + b_n) 
!                              = inv(A) RHS = X
! 
!     INPUT       MEANING
!     -----       ------- 
!     IBPHASE     control on backpropgater correction term, see function hdiffer
!     FREQ        current frequency (Hz) 
!     MDIM        leading dimension
!     MRDOF       holds receiver degree of freedom numbers
!     NDIM        number of components in solution
!     NDOF        number of degrees of freedom
!     NREC        number of receivers 
!     RESID       preconditioned residuals for backpropagation on receiver components
!      
!     OUTPUT      MEANING 
!     ------      -------
!     IRHS_SPARSE holds each RHSs DOF number
!     IRHS_PTR    pointer vector for IRHS_SPARSE 
!     IWARN       warning flag, probably something really wrong 
!     RHS_SPARSE  sparse RHS of residual data to backpropgate  
!
!.... variable declarations
      COMPLEX*8, INTENT(IN) :: RESID(MDIM,*)
      REAL*8, INTENT(IN) :: FREQ
      INTEGER*4, INTENT(IN) :: MRDOF(MDIM,*), MDIM,NDOF,NREC,NDIM, IBPHASE 
      COMPLEX*8, INTENT(OUT) :: RHS_SPARSE(NDIM*NREC)
      INTEGER*4, INTENT(OUT) :: IRHS_SPARSE(NDIM*NREC),IRHS_PTR(2)
!.... local variables
      COMPLEX*8 HFACT, HDIFFER, CZERO
      PARAMETER(CZERO = CMPLX(0.0,0.0))
! 
!----------------------------------------------------------------------------------------!
! 
!.... calculate geometric spreading correction
      IWARN = 0
      HFACT = HDIFFER(IBPHASE,FREQ)
! 
!.... stack on receivers where each receiver acts as a point source 
      IRHS_PTR(1) = 1
      RHS_SPARSE(1:NDIM*NREC) = CZERO
      IZRHS = 0
      DO 2 IREC=1,NREC
         DO 3 I=1,NDIM
! 
!.......... overlay sources using superposition 
            IDOF = MRDOF(I,IREC)
            IF (IDOF <= 0 .OR. IDOF > NDOF) THEN
               IWARN = 1
               WRITE(*,*) 'resid_rhs: Invalid degree of freedom',IDOF,' skipping...'
            ELSE
               IZRHS = IZRHS + 1
               RHS_SPARSE(IZRHS) = RHS_SPARSE(IZRHS) + HFACT*RESID(I,IREC)
               IRHS_SPARSE(IZRHS) = IDOF
            ENDIF
    3    CONTINUE !loop on components in solution 
    2 CONTINUE
      NZRHS = IZRHS
      IRHS_PTR(2) = NZRHS + 1
! 
!.... need to conjugate residuals 
      DO 4 IZRHS=1,NZRHS
         RHS_SPARSE(IZRHS) = CONJG(RHS_SPARSE(IZRHS))
    4 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE RESID_DRHS(MDIM,NDOF,NREC,NDIM, IBPHASE,FREQ,MRDOF, RESID,  RHS,IWARN)
! 
!     This routine calculate the residual populating the RHS for backpropgation.  Here, 
!     we use the principal of superposition to force our RHS into one vector.   
! 
!        x_1 = inv(A) b_1
!        x_2 = inv(A) b_2
!          . 
!          . 
!          .  
!        x_n = inv(A) b_n
! 
!     Then sum to obtain 
!        x_1 + x_2 + ... + x_n = inv(A) b_1 + inv(A) b_2 + ... inv(A) b_n 
!                              = inv(A) (b_1 + b_2 + ... + b_n) 
!                              = inv(A) RHS = X
! 
!     INPUT       MEANING
!     -----       ------- 
!     IBPHASE     control on backpropgater correction term, see function hdiffer
!     FREQ        current frequency (Hz) 
!     MDIM        leading dimension
!     MRDOF       holds receiver degree of freedom numbers
!     NDIM        number of components in solution
!     NDOF        number of degrees of freedom
!     NREC        number of receivers 
!     RESID       preconditioned residuals for backpropagation on receiver components
!      
!     OUTPUT      MEANING 
!     ------      -------
!     IWARN       warning flag, probably something really wrong 
!     RHS         RHS of residual data to backpropgate  
!
!.... variable declarations
      COMPLEX*8, INTENT(IN) :: RESID(MDIM,*)
      REAL*8, INTENT(IN) :: FREQ
      INTEGER*4, INTENT(IN) :: MRDOF(MDIM,*), MDIM,NDOF,NREC,NDIM, IBPHASE
      COMPLEX*8, INTENT(OUT) :: RHS(NDOF)
      INTEGER*4, INTENT(OUT) :: IWARN

!.... local variables
      COMPLEX*8 HFACT, HDIFFER, CZERO
      PARAMETER(CZERO = CMPLX(0.0,0.0))
! 
!----------------------------------------------------------------------------------------!
! 
!.... calculate geometric spreading correction
      IWARN = 0 
      HFACT = HDIFFER(IBPHASE,FREQ)
! 
!.... null out the RHS 
      DO 1 IDOF=1,NDOF
         RHS(IDOF) = CZERO
    1 CONTINUE
!
!.... stack each RHS where each observation is a point source 
      DO 2 IREC=1,NREC
         DO 3 I=1,NDIM
! 
!.......... overlay sources using superposition 
            IDOF = MRDOF(I,IREC)
            IF (IDOF <= 0 .OR. IDOF > NDOF) THEN
               IWARN = 1 
               WRITE(*,*) 'resid_drhs: Invalid degree of freedom',IDOF,' skipping...'
            ELSE
               RHS(IDOF) = RHS(IDOF) + HFACT*RESID(I,IREC)
            ENDIF
    3    CONTINUE !loop on components in solution 
    2 CONTINUE
! 
!.... need to conjugate residuals for solution phase 
      DO 4 IDOF=1,NDOF
         RHS(IDOF) = CONJG(RHS(IDOF))
    4 CONTINUE
      RETURN
      END 

!                                                                                        !
!========================================================================================!
!                                                                                        !
      COMPLEX*8 FUNCTION HDIFFER(IBPHASE,FREQ)
!----------------------------------------------------------------------------------------!
!     Apply a geometric spreading compensation at frequency freq (Hz)                    !
!       ibphase = 0     No shift                                                         !
!       ibphase = 1     Standard half differentiator  sqrt(-iw)                          !
!       ibphase = 2     180 degree phase shifted half differentiator -sqrt(iw)           !
!       ibphase = 3     Conjugate half differentiator sqrt(iw)                           !
!       ibphase = 4     180 degree phase shifted conjugate half differentiator -sqrt(iw) !
!       ibphase = 5     Full differentiator (iw)                                         !
!       ibphase = 6     Hilbert transformer (i)                                          !
!----------------------------------------------------------------------------------------!
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8 FREQ
      REAL*4 OMEGA
      INTEGER*4 IBPHASE
      REAL*8 PI/3.1415926535897932384626434D0/
!
!----------------------------------------------------------------------------------------!
!
!.... calculate the half differentiator
      OMEGA = REAL(2.D0*PI*FREQ)
      HDIFFER = CMPLX(1.0,0.0)
      IF (IBPHASE == 1 .OR. IBPHASE == 2) THEN
         HDIFFER = SQRT(CMPLX(0.0,-OMEGA))
         IF (IBPHASE == 2) HDIFFER =-HDIFFER
      ELSEIF (IBPHASE == 3 .OR. IBPHASE == 4) THEN
         HDIFFER = SQRT(CMPLX(0.0, OMEGA))
         IF (IBPHASE == 4) HDIFFER =-HDIFFER
      ELSEIF (IBPHASE == 5) THEN
         HDIFFER = CMPLX(0.0,OMEGA)
      ELSEIF (IBPHASE == 6) THEN
         HDIFFER = CMPLX(0.0,1.0)
      ELSE
         HDIFFER = CMPLX(1.0,0.0)
      ENDIF
      RETURN
      END
