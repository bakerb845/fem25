      SUBROUTINE GET_OE_SB(LSURF, MDIM,MFREQ,MREC, INV,FRQ,SRC,RCV,  &
                           WGT,EST,OBS,IERR) 
!
!     Get observations or estimates corresponding to a surface or body wave 
!     frequency at modeling frequencies.  Useful for STF and RRF updates
!
!     INPUT      MEANING
!     -----      -------
!     LSURF      True -> get surface wave data / estimates
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!
!     OUTPUT     MEANING
!     ------     ------- 
!     EST        estimates corresponding to only surface or body wave modeling 
!                frequencies
!     IERR       error flag
!     OBS        observations corresponding to only surface or body wave modeling
!                frequencies
!
!.... variable declarations 
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE(SRC_INFO)  SRC  
      TYPE(RECV_INFO) RCV
      TYPE(FRQ_INFO)  FRQ  
      TYPE(INV_INFO)  INV  
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MREC 
      LOGICAL*4, INTENT(IN) :: LSURF
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MFREQ,MREC,*), OBS(MDIM,MFREQ,MREC,*)
      REAL*4, INTENT(OUT) :: WGT(MDIM,MFREQ,MREC,*)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      INTEGER*4 IFREQ, IFREQ_LOC, ISRC, ISRC_LOC, IREC, I
      LOGICAL*4 LBODY 
!
!----------------------------------------------------------------------------------------!
!
!.... determine if we want body or surface waves
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
      IFREQ_LOC = 0
      DO 1 IFREQ=1,frq%NFREQ
         IF (LSURF) THEN 
            IF (frq%CFTYPE(IFREQ) == 'B') GOTO 100
         ENDIF
         IF (LBODY) THEN 
            IF (frq%CFTYPE(IFREQ) == 'S') GOTO 100
         ENDIF 
         IFREQ_LOC = IFREQ_LOC + 1
         ISRC_LOC = 0
         DO 2 ISRC=1,src%NSRC
            IF (frq%CFTYPE(IFREQ) == 'S' .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
                frq%CFTYPE(IFREQ) == 'B' .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN 
               ISRC_LOC = ISRC_LOC + 1
               DO 3 IREC=1,rcv%NREC
                  DO 4 I=1,NDIM
                     EST(I,IFREQ_LOC,IREC,ISRC_LOC) =   inv%EST(I,IFREQ,IREC,ISRC)
                     OBS(I,IFREQ_LOC,IREC,ISRC_LOC) =   inv%OBS(I,IFREQ,IREC,ISRC)
                     WGT(I,IFREQ_LOC,IREC,ISRC_LOC) = inv%WGHTS(I,IFREQ,IREC,ISRC)
    4             CONTINUE !loop on components
    3          CONTINUE !loop on receiver
            ENDIF !end check on source frequency match
    2    CONTINUE !loop on sources
!
!....... error checks
         IF (LSURF) THEN
            IF (ISRC_LOC /= src%NSRC_SRF) THEN
               WRITE(*,*) 'get_oe_sb: Error isrc_loc /= src%nsrc_srf!'
               IERR = 1
               RETURN
            ENDIF
         ELSE  !body wave
            IF (ISRC_LOC /= src%NSRC_BDY) THEN
               WRITE(*,*) 'get_oe_sb: Error isrc_loc /= src%nsrc_bdy!'
               IERR = 1
               RETURN
            ENDIF 
         ENDIF
  100    CONTINUE
    1 CONTINUE
!
!.... error checks
      IF (LSURF) THEN
         IF (IFREQ_LOC /= frq%NFREQ_SRF) THEN
            WRITE(*,*) 'get_oe_sb: Error ifreq_loc /= frq%nfreq_srf'
            IERR = 2
         ENDIF
      ELSE
         IF (IFREQ_LOC /= frq%NFREQ_BDY) THEN
            WRITE(*,*) 'get_oe_sb: Error ifreq_loc /= frq%nfreq_bdy'
            IERR = 2
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_EST_SB(LSURF, MDIM,MFREQ,MREC, INV,FRQ,SRC,RCV, EST,IERR) 
!
!     Get estimates corresponding to a surface or body wave 
!     frequency at modeling frequencies.  Useful for STF and RRF updates and
!     file IO
!
!     INPUT      MEANING
!     -----      -------
!     LSURF      True -> get surface wave data / estimates
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!
!     OUTPUT     MEANING
!     ------     ------- 
!     EST        estimates corresponding to only surface or body wave modeling 
!                frequencies
!     IERR       error flag
!
!.... variable declarations 
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE(SRC_INFO)  SRC
      TYPE(RECV_INFO) RCV
      TYPE(FRQ_INFO)  FRQ
      TYPE(INV_INFO)  INV
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MREC
      LOGICAL*4, INTENT(IN) :: LSURF
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MFREQ,MREC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 IFREQ, IFREQ_LOC, ISRC, ISRC_LOC, IREC, I
      LOGICAL*4 LBODY
!
!----------------------------------------------------------------------------------------!
!
!.... determine if we want body or surface waves
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
      IFREQ_LOC = 0
      DO 1 IFREQ=1,frq%NFREQ
         IF (LSURF) THEN
            IF (frq%CFTYPE(IFREQ) == 'B') GOTO 100
         ENDIF
         IF (LBODY) THEN
            IF (frq%CFTYPE(IFREQ) == 'S') GOTO 100
         ENDIF
         IFREQ_LOC = IFREQ_LOC + 1
         ISRC_LOC = 0
         DO 2 ISRC=1,src%NSRC
            IF (frq%CFTYPE(IFREQ) == 'S' .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
                frq%CFTYPE(IFREQ) == 'B' .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
               ISRC_LOC = ISRC_LOC + 1
               DO 3 IREC=1,rcv%NREC
                  DO 4 I=1,NDIM
                     EST(I,IFREQ_LOC,IREC,ISRC_LOC) =   inv%EST(I,IFREQ,IREC,ISRC)
    4             CONTINUE !loop on components
    3          CONTINUE !loop on receiver
            ENDIF !end check on source frequency match
    2    CONTINUE !loop on sources
!
!....... error checks
         IF (LSURF) THEN
            IF (ISRC_LOC /= src%NSRC_SRF) THEN
               WRITE(*,*) 'get_est_sb: Error isrc_loc /= src%nsrc_srf!'
               IERR = 1
               RETURN
            ENDIF
         ELSE  !body wave
            IF (ISRC_LOC /= src%NSRC_BDY) THEN
               WRITE(*,*) 'get_est_sb: Error isrc_loc /= src%nsrc_bdy!'
               IERR = 1
               RETURN
            ENDIF
         ENDIF
  100    CONTINUE
    1 CONTINUE
!
!.... error checks
      IF (LSURF) THEN 
         IF (IFREQ_LOC /= frq%NFREQ_SRF) THEN 
            WRITE(*,*) 'get_est_sb: Error ifreq_loc /= frq%nfreq_srf'
            IERR = 2
         ENDIF
      ELSE 
         IF (IFREQ_LOC /= frq%NFREQ_BDY) THEN 
            WRITE(*,*) 'get_est_sb: Error ifreq_loc /= frq%nfreq_bdy'
            IERR = 2
         ENDIF
      ENDIF
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_TWIN(MREC,MSRC, NREC, LSURF,NSRC_REF, WIN,SRC, TWIN,IERR) 
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (SRC_INFO) SRC
      TYPE (WIN_INFO) WIN 
      INTEGER*4, INTENT(IN) :: MREC, MSRC, NREC, NSRC_REF 
      LOGICAL*4, INTENT(IN) :: LSURF 
      REAL*4, INTENT(OUT) :: TWIN(MREC,MSRC,*)  
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 ISRC, JSRC, IREC
      LOGICAL*4 LBODY 
!
!----------------------------------------------------------------------------------------!
!
!.... determine if we want body or surface waves
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on sources and extract
      JSRC = 0 
      DO 1 ISRC=1,src%NSRC
         IF (LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR. &
             LBODY .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
            JSRC = JSRC + 1
            DO 2 IREC=1,NREC
               TWIN(IREC,JSRC,1) = win%TWIN(IREC,ISRC,1)
               TWIN(IREC,JSRC,2) = win%TWIN(IREC,ISRC,2)
    2       CONTINUE
         ENDIF
    1 CONTINUE !loop on sources
      IF (JSRC /= NSRC_REF) THEN
         IERR = 1
         WRITE(*,*) 'get_twin: Error jsrc /= nsrc_ref!'
         RETURN
      ENDIF 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SET_EST_SB(LSURF, MDIM,MFREQ,MREC, EST, INV,FRQ,SRC,RCV, IERR) 
!
!     Puts the local estimates back onto the global  estimates
!
!     INPUT      MEANING
!     -----      -------
!     EST        local estaimtes for surface or body wave modeling
!     LSURF      True -> get surface wave data / estimates
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!
!.... variable declarations 
      IMPLICIT NONE 
      INCLUDE 'fwd_struc.h'
      TYPE(SRC_INFO)  SRC  
      TYPE(RECV_INFO) RCV
      TYPE(FRQ_INFO)  FRQ  
      TYPE(INV_INFO)  INV  
      COMPLEX*8, INTENT(IN) :: EST(MDIM,MFREQ,MREC,*)
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MREC 
      LOGICAL*4, INTENT(IN) :: LSURF
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      INTEGER*4 IFREQ, IFREQ_LOC, ISRC, ISRC_LOC, IREC, I
      LOGICAL*4 LBODY 
!
!----------------------------------------------------------------------------------------!
!
!.... determine if we want body or surface waves
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
      IFREQ_LOC = 0
      DO 1 IFREQ=1,frq%NFREQ
         IF (LSURF) THEN
            IF (frq%CFTYPE(IFREQ) == 'B') GOTO 100
         ENDIF
         IF (LBODY) THEN
            IF (frq%CFTYPE(IFREQ) == 'S') GOTO 100
         ENDIF
         IFREQ_LOC = IFREQ_LOC + 1
         ISRC_LOC = 0
         DO 2 ISRC=1,src%NSRC
            IF (frq%CFTYPE(IFREQ) == 'S' .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
                frq%CFTYPE(IFREQ) == 'B' .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
               ISRC_LOC = ISRC_LOC + 1
               DO 3 IREC=1,rcv%NREC
                  DO 4 I=1,NDIM
                     inv%EST(I,IFREQ,IREC,ISRC) = EST(I,IFREQ_LOC,IREC,ISRC_LOC) 
    4             CONTINUE !loop on components
    3          CONTINUE !loop on receiver
            ENDIF !end check on source frequency match
    2    CONTINUE !loop on sources
!
!....... error checks
         IF (LSURF) THEN
            IF (ISRC_LOC /= src%NSRC_SRF) THEN
               WRITE(*,*) 'get_oe_sb: Error isrc_loc /= src%nsrc_srf!'
               IERR = 1
               RETURN
            ENDIF
         ELSE  !body wave
            IF (ISRC_LOC /= src%NSRC_BDY) THEN
               WRITE(*,*) 'get_oe_sb: Error isrc_loc /= src%nsrc_bdy!'
               IERR = 1
               RETURN
            ENDIF
         ENDIF
  100    CONTINUE
    1 CONTINUE
!
!.... error checks
      IF (LSURF) THEN
         IF (IFREQ_LOC /= frq%NFREQ_SRF) THEN
            WRITE(*,*) 'get_oe_sb: Error ifreq_loc /= frq%nfreq_srf'
            IERR = 2
         ENDIF
      ELSE
         IF (IFREQ_LOC /= frq%NFREQ_BDY) THEN
            WRITE(*,*) 'get_oe_sb: Error ifreq_loc /= frq%nfreq_bdy'
            IERR = 2
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_OEINV_SB(LSURF, MDIM,MFREQ,MREC, INV,FRQ,SRC,RCV,  &
                              WGT,EST,OBS,IERR) 
!
!     Get observations or estimates corresponding to a surface or body wave 
!     frequency at inversion frequencies only.  Useful for inversion
!
!     INPUT      MEANING
!     -----      -------
!     LSURF      True -> get surface wave data / estimates
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!
!     OUTPUT     MEANING
!     ------     ------- 
!     EST        estimates corresponding to only surface or body wave inversion 
!                frequencies
!     IERR       error flag
!     OBS        observations corresponding to only surface or body wave inversion
!                frequencies
!
!.... variable declarations 
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE(SRC_INFO)  SRC  
      TYPE(RECV_INFO) RCV 
      TYPE(FRQ_INFO)  FRQ  
      TYPE(INV_INFO)  INV  
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MREC 
      LOGICAL*4, INTENT(IN) :: LSURF
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MFREQ,MREC,*), OBS(MDIM,MFREQ,MREC,*)
      REAL*4, INTENT(OUT) :: WGT(MDIM,MFREQ,MREC,*)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      INTEGER*4 IFREQ, IFREQ_LOC, ISRC, ISRC_LOC, IREC, I
      LOGICAL*4 LBODY 
!
!----------------------------------------------------------------------------------------!
!
!.... determine if we want body or surface waves
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
      IFREQ_LOC = 0
      DO 1 IFREQ=1,frq%NFREQ
         IF (LSURF) THEN
            IF (frq%CFTYPE(IFREQ) == 'B') GOTO 100
            IF (.NOT.frq%LINVF(IFREQ)) GOTO 100 
         ENDIF
         IF (LBODY) THEN
            IF (frq%CFTYPE(IFREQ) == 'S') GOTO 100
            IF (.NOT.frq%LINVF(IFREQ)) GOTO 100
         ENDIF
         IFREQ_LOC = IFREQ_LOC + 1
         ISRC_LOC = 0
         DO 2 ISRC=1,src%NSRC
            IF (frq%CFTYPE(IFREQ) == 'S' .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
                frq%CFTYPE(IFREQ) == 'B' .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
               ISRC_LOC = ISRC_LOC + 1
               DO 3 IREC=1,rcv%NREC
                  DO 4 I=1,NDIM
                     EST(I,IFREQ_LOC,IREC,ISRC_LOC) =   inv%EST(I,IFREQ,IREC,ISRC)
                     OBS(I,IFREQ_LOC,IREC,ISRC_LOC) =   inv%OBS(I,IFREQ,IREC,ISRC)
                     WGT(I,IFREQ_LOC,IREC,ISRC_LOC) = inv%WGHTS(I,IFREQ,IREC,ISRC)
    4             CONTINUE !loop on components
    3          CONTINUE !loop on receiver
            ENDIF !end check on source frequency match
    2    CONTINUE !loop on sources
!
!....... error checks
         IF (LSURF) THEN
            IF (ISRC_LOC /= src%NSRC_SRF) THEN
               WRITE(*,*) 'get_oeinv_sb: Error isrc_loc /= src%nsrc_srf!'
               IERR = 1
               RETURN
            ENDIF
         ELSE  !body wave
            IF (ISRC_LOC /= src%NSRC_BDY) THEN
               WRITE(*,*) 'get_oeinv_sb: Error isrc_loc /= src%nsrc_bdy!'
               IERR = 1
               RETURN
            ENDIF
         ENDIF
  100    CONTINUE
    1 CONTINUE
!
!.... error checks
      IF (LSURF) THEN
         IF (IFREQ_LOC /= frq%NFREQ_SRF_INV) THEN
            WRITE(*,*) 'get_oeinv_sb: Error ifreq_loc /= frq%nfreq_srf_inv'
            IERR = 2
         ENDIF
      ELSE
         IF (IFREQ_LOC /= frq%NFREQ_BDY_INV) THEN
            WRITE(*,*) 'get_oeinv_sb: Error ifreq_loc /= frq%nfreq_bdy_inv'
            IERR = 2
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE MULJAC_RESID_MPI(MYID,MASTER,MYCOMM, NPPGRP,             &
                                  LSURF,LSTF, AZMOD, INV,FRQ,SRC,RCV, IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO)  INV 
      TYPE (FRQ_INFO)  FRQ
      TYPE (SRC_INFO)  SRC
      TYPE (RECV_INFO) RCV
      REAL*8, INTENT(IN) :: AZMOD
      INTEGER*4, INTENT(IN) :: MYID, MASTER, MYCOMM, NPPGRP
      LOGICAL*4, INTENT(IN) :: LSURF
      LOGICAL*4, INTENT(INOUT) :: LSTF 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*8, ALLOCATABLE :: OBS(:,:,:,:), EST(:,:,:,:), SRCWRK(:,:)
      REAL*8, ALLOCATABLE :: COVD8(:,:), PYLOC(:,:), FREQ_LOC(:)
      REAL*4, ALLOCATABLE :: WGT(:,:,:,:), JMATL(:,:), W1(:,:), COVD(:,:), WORK(:), &
                             RESIDV(:), BUFF(:), GRAD(:)
      INTEGER*4, ALLOCATABLE :: IPERM_REC(:)
      CHARACTER(80) FLNAME
      CHARACTER(12) CFREQ
      CHARACTER(5)  CSRC, CID
      COMPLEX*8 CTEMP
      REAL*8 TEMP
      REAL*4 ALPHA, OBJBUF, DTMAX, OBJ, SDOT
      INTEGER*4 NOBS2, NOBS, LDJAC, N, NFREQ, NSRC, IFREQ, IFREQL, ISRC, IREC, JSRC,    &
                LDC, M, IPPGRP, I, J, JREC, I1, I2, ISZERO, NREC, ML, MYRANK,  &
                IROW, JCOL, LROW, KCOL, INDX, JNDX, JFREQ, NPROCS, MPIERR
      LOGICAL*4 LBODY
      INTEGER*4 NOBSLST
      REAL*8, PARAMETER :: TWOPI = 6.2831853071795862
      LOGICAL*4, PARAMETER :: LSHIFT = .TRUE.
!
!----------------------------------------------------------------------------------------!
!
!.... sort ou MPI information
      CALL MPI_COMM_RANK(MYCOMM,MYRANK,MPIERR)
      CALL MPI_COMM_SIZE(MYCOMM,NPROCS,MPIERR)
!.... broadcast the residuals and calculate the gradient
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... get observations and estimates for surface or body wave
      NREC = rcv%NREC
      NOBS2 = 2*NDIM*NREC 
      LDJAC = NOBS2
      N = inv%NA35 
      IF (LSURF) THEN
         NFREQ = frq%NFREQ_SRF_INV
         NSRC  = src%NSRC_SRF
         ALLOCATE(OBS(NDIM,NFREQ,NREC,NSRC))
         ALLOCATE(EST(NDIM,NFREQ,NREC,NSRC))
         ALLOCATE(WGT(NDIM,NFREQ,NREC,NSRC))
         OBS(:,:,:,:) = CMPLX(0.0,0.0)
         EST(:,:,:,:) = CMPLX(0.0,0.0)
         WGT(:,:,:,:) = 0.0
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'muljac_resid_mpi: Caclulating surface wave gradient...'
            CALL GET_OEINV_SB(.TRUE., NDIM,NFREQ,NREC, INV,FRQ,SRC,RCV,  &
                              WGT,EST,OBS,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'muljac_resid_mpi: Error calling get_oeinv_sb 1!'
               RETURN
            ENDIF
         ENDIF 
         DTMAX = inv%DTMAX_SRF
      ELSE
         NFREQ = frq%NFREQ_BDY_INV
         NSRC  = src%NSRC_BDY
         ALLOCATE(OBS(NDIM,NFREQ,NREC,NSRC))
         ALLOCATE(EST(NDIM,NFREQ,NREC,NSRC))
         ALLOCATE(WGT(NDIM,NFREQ,NREC,NSRC))
         OBS(:,:,:,:) = 0.0
         EST(:,:,:,:) = 0.0
         WGT(:,:,:,:) = 0.0
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'muljac_resid_mpi: Calculating body wave gradient...'
            CALL GET_OEINV_SB(.FALSE., NDIM,NFREQ,NREC, INV,FRQ,SRC,RCV, &
                              WGT,EST,OBS,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'muljac_resid_mpi: Error calling get_oeinv_sb 2!'
               RETURN
            ENDIF
         ENDIF
         DTMAX = inv%DTMAX_BDY
      ENDIF
      M = 2*NDIM*rcv%NREC
      LDC = M
      ALLOCATE(COVD8(LDC,LDC))
      COVD8(:,:) = 0.D0
!     IF (LSURF .AND. inv%LCOVD_SRF .OR. LBODY .AND. inv%LCOVD_BDY) THEN
       IF (LSURF .AND. .NOT.inv%LCOVD_SRF .OR. &
           LBODY .AND. .NOT.inv%LCOVD_BDY) THEN 
         COVD8(:,:) = 0.D0
         DO I=1,M
            COVD8(I,I) = 1.D0 
         ENDDO
      ENDIF
      ALLOCATE(FREQ_LOC(NFREQ)) 
      CALL GET_FREQ_INV(frq%NFREQ, NFREQ, LSURF, frq%CFTYPE,frq%LINVF,frq%FREQ,  &
                        FREQ_LOC,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'muljac_resid_mpi: Error calling get_freq_inv!'
         RETURN
      ENDIF
      allocate(pyloc(nfreq,nsrc))
      jfreq = 0
      do ifreq=1,frq%nfreq
         if (     lsurf .and. frq%cftype(ifreq) == 'S' .and. frq%linvf(ifreq) .or. &
             .not.lsurf .and. frq%cftype(ifreq) == 'B' .and. frq%linvf(ifreq)) then 
            jfreq = jfreq + 1
            jsrc = 0
            do isrc=1,src%nsrc
               if (     lsurf .and. src%srctyp(isrc)(1:1) == 'S' .or. &
                   .not.lsurf .and. src%srctyp(isrc)(1:1) == 'P') then
                  jsrc = jsrc + 1 
                  pyloc(jfreq,jsrc) = src%pytab(ifreq,isrc)  
               endif
            enddo
         endif
      enddo
!
!.... broadcast
      IF (MYID == MASTER) WRITE(*,*) 'muljac_resid_mpi: Broadcasting...'
      CALL MPI_BCAST(LSTF,1,MPI_LOGICAL,MASTER,MYCOMM,MPIERR)
      IF (LSTF) THEN
         ALLOCATE(SRCWRK(NFREQ,NSRC))
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'muljac_resid_mpi: Fetching STF...'
            IF (LSURF) THEN
               CALL GET_SRC_INV(frq%NFREQ,NFREQ, frq%NFREQ,src%NSRC,  &
                                frq%NFREQ_SRF_INV,src%NSRC_SRF, .TRUE., &
                                frq%CFTYPE, src%SRCTYP, frq%LINVF, src%SOURCE,  &
                                SRCWRK, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'muljac_resid_mpi: Error calling get_src_inv 1'
                  RETURN
               ENDIF
            ELSE !body wave
               CALL GET_SRC_INV(frq%NFREQ,NFREQ, frq%NFREQ,src%NSRC,  &
                                frq%NFREQ_BDY_INV,src%NSRC_BDY, .FALSE., &
                                frq%CFTYPE, src%SRCTYP, frq%LINVF, src%SOURCE,  &
                                SRCWRK, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'muljac_resid_mpi: Error calling get_src_inv 2'
                  RETURN
               ENDIF
            ENDIF ! end check on body or surface wave
         ENDIF !end check on myid
      ENDIF !end check on if we need a source time function
      DO 1 IFREQ=1,NFREQ
         DO 2 ISRC=1,NSRC
            DO 3 I=1,NDIM
               CALL MPI_BCAST(WGT(I,IFREQ,1:NREC,ISRC),NREC,MPI_REAL,    &
                              MASTER,MYCOMM,MPIERR)
               CALL MPI_BCAST(OBS(I,IFREQ,1:NREC,ISRC),NREC,MPI_COMPLEX, &
                              MASTER,MYCOMM,MPIERR)
               CALL MPI_BCAST(EST(I,IFREQ,1:NREC,ISRC),NREC,MPI_COMPLEX, &
                              MASTER,MYCOMM,MPIERR)
    3       CONTINUE
            IF (LSTF) THEN
               CALL MPI_BCAST(SRCWRK(IFREQ,ISRC),1,MPI_COMPLEX, &
                              MASTER,MYCOMM,MPIERR)
            ENDIF
    2    CONTINUE
    1 CONTINUE 
!
!.... set space
      ALLOCATE(W1(NDIM,rcv%NREC))
      ALLOCATE(COVD(LDC,LDC))
      ALLOCATE(RESIDV(NOBS2)) 
      ALLOCATE(WORK(NOBS2))
      ALLOCATE(JMATL(LDJAC,N))  
      ALLOCATE(BUFF(N))
      ALLOCATE(IPERM_REC(NDIM*NREC))
      BUFF(1:N) = 0.0 
      OBJBUF = 0.0
      JMATL(:,:) = 0.0
!
!.... loop on frequencies
      ISZERO = 0
!     if (myid > 0) goto 1110
      ISZERO = 1
      DO 11 IFREQL=1,NFREQ,NPROCS
         IFREQ = IFREQL + MYRANK
         IF (IFREQ > NFREQ) GOTO 101     !out of inversion frequencies
         CFREQ(1:12) = ' '
         WRITE(CFREQ,'(F12.5)') FREQ_LOC(IFREQ) !frq%FREQ(IFREQ)
         CFREQ = ADJUSTL(CFREQ)
!
!....... loop on sources
         JSRC = 0
         DO 12 ISRC=1,src%NSRC  
            IF (LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'P') GOTO 120  
            IF (LBODY .AND. src%SRCTYP(ISRC)(1:1) == 'S') GOTO 120
            CSRC(1:5) = ' '
            WRITE(CSRC,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC)
            JSRC = JSRC + 1 
            WORK(1:NOBS2) = 0.0
            RESIDV(1:NOBS2) = 0.0 
            COVD(:,:) = 0.0
            NOBS = NOBSLST(NDIM, NREC,NDIM,  OBS(1:NDIM,IFREQ,1:NREC,JSRC)) 
            IF (NOBS > 0) THEN
               CALL PERM_OBS(NDIM,NDIM,NREC, OBS(1:NDIM,IFREQ,1:NREC,JSRC), &
                             IPERM_REC,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'muljac_resid_mpi: Error calling perm_obs!'
                  RETURN
               ENDIF 
               IF (LSURF.AND.inv%LCOVD_SRF .OR. LBODY.AND.inv%LCOVD_BDY) THEN 
                  CALL LAPLACE_1Dv2(LDC,NDIM,NREC,inv%LDWGHT, inv%DX,  &
                                    OBS(1:NDIM,IFREQ,1:NREC,JSRC), COVD8,IERR)
                  IF (IERR /= 0) THEN 
                     WRITE(*,*) 'muljac_resid_mpi: Error calling laplace_1d!'
                     COVD8(:,:) = 0.D0
                     DO 9 I=1,M
                        COVD8(I,I) = 1.D0 
    9                CONTINUE
                  ENDIF
               ENDIF

               DO 13 IREC=1,NREC
                  DO 14 I=1,NDIM
                     IF (WGT(I,IFREQ,IREC,JSRC) > 0.0) THEN
                        W1(I,IREC) = 1.0
                     ELSE
                        W1(I,IREC) = 0.0
                     ENDIF
   14             CONTINUE
   13          CONTINUE
               CALL OE2RESV(NDIM, NOBS2, NDIM,NREC, DTMAX,                   &
                            inv%IRESTP,0, AZMOD,FREQ_LOC(IFREQ),     &
                            IPERM_REC, W1, &
                            EST(1:NDIM,IFREQ,1:NREC,JSRC),  &
                            OBS(1:NDIM,IFREQ,1:NREC,JSRC), RESIDV,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'muljac_resid_mpi: Error calling oe2resv!'
                  RETURN
               ENDIF
               IF (LSHIFT) THEN
                  DO IREC=1,rcv%NREC
                     DO I=1,NDIM 
                        INDX =             IPERM_REC((IREC - 1)*NDIM + I) 
                        JNDX = NDIM*NREC + IPERM_REC((IREC - 1)*NDIM + I)
                        TEMP  = TWOPI*FREQ_LOC(IFREQ)*PYLOC(IFREQ,JSRC)*rcv%YREC(IREC)
                        CTEMP = CMPLX(RESIDV(INDX),RESIDV(JNDX))
                        !phase shift back into plane -i w py y
                        CTEMP = CEXP(CMPLX(0.0,-SNGL(TEMP)))*CTEMP
                        RESIDV(INDX) = REAL(CTEMP)
                        RESIDV(JNDX) = IMAG(CTEMP)
                     ENDDO
                  ENDDO
               ENDIF
!
!............. scale data covariance matrix by weights
               I = 1
               IREC = 0
               COVD(:,:) = 0.0
               DO 15 IROW=1,M
                  IREC = IREC + 1
                  IF (IREC > NREC) THEN
                     IREC = 1
                     I = I + 1
                  ENDIF
                  IF (IROW == NDIM*rcv%NREC + 1) I = 1
                  IF (IROW > NDIM*rcv%NREC) THEN
                     LROW = NDIM*NREC + IPERM_REC( (IREC - 1)*NDIM + I)
                  ELSE
                     LROW = IPERM_REC( (IREC - 1)*NDIM + I)
                  ENDIF
                  J = 1
                  JREC = 0
                  DO 16 JCOL=1,M
                     JREC = JREC + 1
                     IF (JREC > NREC) THEN
                        JREC = 1
                        J = J + 1
                     ENDIF
                     IF (JCOL == NDIM*rcv%NREC + 1) J = 1
                     IF (JCOL > NDIM*rcv%NREC) THEN
                         KCOL = NDIM*rcv%NREC + IPERM_REC( (JREC - 1)*NDIM + J)
                     ELSE
                         KCOL = IPERM_REC( (JREC - 1)*NDIM + J)
                     ENDIF
                     COVD(LROW,KCOL) = SQRT(WGT(I,IFREQ,IREC,JSRC)) &
                                      *SNGL(COVD8(LROW,KCOL)) &
                                      *SQRT(WGT(J,IFREQ,JREC,JSRC))
   16             CONTINUE 
   15          CONTINUE
!
!............. calculate Sigma^{-1} C_D^{-1} Sigma^{-1} (d - Au) 
!              WORK(:) = 0.0
!              CALL SCOPY(M,RESIDV,1,WORK,1)
!              RESIDV(:) = 0.0
!              CALL SGEMV('N',M,M, 1.0, COVD,LDC, WORK,1, 0.0,RESIDV,1)
!
!............. loop on process groups 
               DO 17 IPPGRP=1,NPPGRP
                  CID(1:5) = ' '
                  WRITE(CID,'(I5)') IPPGRP - 1
                  CID = ADJUSTL(CID)
                  FLNAME(1:80) = ' '
                  IF (LSURF) THEN
                     FLNAME = './scratch/jac_surf_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'//  &
                              TRIM(CSRC)//'.dat'
                  ELSE
                     FLNAME = './scratch/jac_body_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'//  &
                             TRIM(CSRC)//'.dat'
                  ENDIF
                  FLNAME = ADJUSTL(FLNAME)
                  CALL LOAD_FJAC(FLNAME,LDJAC, NDIM,NREC, IPERM_REC, JMATL,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'muljac_resid_mpi: Unable to load Jacobian:',  &
                                TRIM(FLNAME), MYID
                     RETURN
                  ENDIF
   17          CONTINUE !loop on process groups 
            ELSE !no observations -> no jacobian -> break ahead
               GOTO 130
            ENDIF  !end check on nobs > 0
!
!.......... convolve the STF into the Jacobians?
            IF (LSTF) CALL SCALE_JLOC(LDJAC,NREC*NDIM,inv%NA35, &
                                      SRCWRK(IFREQ,JSRC),JMATL)
!
!.......... calculate Sigma^{-1} C_D^{-1} Sigma^{-1} (d - Au) 
            WORK(:) = 0.0
            CALL SCOPY(M,RESIDV,1,WORK,1)
            RESIDV(:) = 0.0
            CALL SGEMV('N',M,M, 1.0, COVD,LDC, WORK,1, 0.0,RESIDV,1)
!           write(48+myid,*) residv
!
!.......... fidelity check on objective function; ((d - u)/s)^2
!           I1 = 0
!           I2 = NDIM*NREC 
!           IOBS = 0
!           DO 20 IREC=1,NREC
!              DO 21 I=1,NDIM 
!                 IOBS = IOBS + 1
!                 I1 = I1 + 1
!                 I2 = I2 + 1
!                 OBJBUF = OBJBUF + RESIDV(I1)**2 + RESIDV(I2)**2
!  21          CONTINUE
!  20       CONTINUE
            OBJBUF = OBJBUF + SDOT(M,RESIDV,1,RESIDV,1)
            !write(48+myid,*) residv
!
!.......... update grad = grad - Re[ adj(J) delta d] =-Re(J^T) Re(d) - Im(J^T) Im(d)
            I1 = 1
            I2 = NDIM*NREC
            ML = I2 - I1 + 1
            ALPHA =-1.0
            CALL SGEMV('T',M,N,ALPHA,JMATL,M, RESIDV,1, 1.0,BUFF,1)
!           CALL SGEMV('T',ML,N,ALPHA,JMATL(I1:I2,1:N),ML,RESIDV(I1:I2),1,1.0,BUFF,1) 
!           I1 = NDIM*NREC + 1
!           I2 = 2*NDIM*NREC
!           ML = I2 - I1 + 1
!           CALL SGEMV('T',ML,N,ALPHA,JMATL(I1:I2,1:N),ML,RESIDV(I1:I2),1,1.0,BUFF,1)
            ISZERO = 0
  130       CONTINUE !no observations
  120       CONTINUE !not the correcto source
   12    CONTINUE !loop on sources
  101    CONTINUE !group out of frequencies
   11 CONTINUE 
!1110 continue
      IERR = 0
      CALL MPI_REDUCE(ISZERO,IERR,1,MPI_INTEGER,MPI_MIN, MASTER,MYCOMM,MPIERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'muljac_resid_MPI: Error the gradient is zero!'
         RETURN
      ENDIF
      CALL MPI_REDUCE(OBJBUF,OBJ,1,MPI_REAL,MPI_SUM, MASTER,MYCOMM,MPIERR)
      IF (MYID == MASTER) THEN
         IF (LSURF) THEN 
            WRITE(*,*) 'muljac_resid_MPI: Surface wave L2 objective function:',OBJ/2.0
         ELSE
            WRITE(*,*) 'muljac_resid_MPI: Body wave L2 objective function:',OBJ/2.0
         ENDIF
      ENDIF
!
!.... reduce onto buffer
      IF (MYID == MASTER) THEN
         ALLOCATE(GRAD(N))
         GRAD(1:N) = 0.0
      ELSE
         ALLOCATE(GRAD(1))
         GRAD(1) = 0.0
      ENDIF
      CALL MPI_REDUCE(BUFF,GRAD,N,MPI_REAL,MPI_SUM, MASTER,MYCOMM,MPIERR)
      IF (MYID == MASTER) THEN
         IF (LSURF) THEN
            CALL SCOPY(N,GRAD,1,inv%GRAD_SRF,1)
         ELSE
            CALL SCOPY(N,GRAD,1,inv%GRAD_BDY,1)
         ENDIF
      ENDIF
!
!.... free space 
      IF (ALLOCATED(FREQ_LOC))  DEALLOCATE(FREQ_LOC) 
      IF (ALLOCATED(RESIDV))    DEALLOCATE(RESIDV)
      IF (ALLOCATED(IPERM_REC)) DEALLOCATE(IPERM_REC)
      IF (ALLOCATED(WGT))       DEALLOCATE(WGT) 
      IF (ALLOCATED(OBS))       DEALLOCATE(OBS)
      IF (ALLOCATED(EST))       DEALLOCATE(EST)
      IF (ALLOCATED(SRCWRK))    DEALLOCATE(SRCWRK)
      IF (ALLOCATED(BUFF))      DEALLOCATE(BUFF) 
      IF (ALLOCATED(GRAD))      DEALLOCATE(GRAD) 
      IF (ALLOCATED(WORK))      DEALLOCATE(WORK)
      IF (ALLOCATED(COVD8))     DEALLOCATE(COVD8)
      IF (ALLOCATED(COVD))      DEALLOCATE(COVD)
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE OE2RESV(MDIM, NOBS2, NDIM,NREC, DTMAX, &
                         IRESTP,IBPHASE,  AZMOD,FREQ,  IPERM_REC, &
                         WGHTS, EST,OBS, RESIDV,IERR)
!
!     Converts the observations and estimates in matrix form to a vector for 
!     matrix vector multiplication, where the Jacobian is stored as 
!
!       J = [ Re{J} ]
!           [ Im{J} ]
!
!     where J is dimension [n_{obs} x n_{inv}].  Furthermore, each row of J is 
!     packed { reciever 1 
!                u,v,w
!              reciever 2
!                u,v,w
!              .
!              .
!              .
!              recevier n_{rec}
!                u,v,w}
!
!     INPUT      MEANING
!     -----      -------
!     AZMOD      model azimuth (deg)
!     DTMAX      max residual (s)
!     EST        complex valued estimates (U,V,W)
!     FREQ       current frequency (Hz)
!     IBPHASE    used if applying half differentiator, pry 0 -> off
!     IFREQ      current frequency number
!     IRESTP     residual type (1) phase, (2) amp, (3) both
!     ISRC       current source number
!     LUNWRAP    True -> need to put obs, est back as complex number
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!     NDIM       number of components (3)
!     NOBS2      number of observations [2*ndim*nrec]
!     NREC       number of recievers
!     OBS        complex valued observations (N,E,Z) 
!     WGHTS      weights on data
!     
!     OUTPUT     MEANING
!     ------     ------ 
!     IERR       error flag
!     RESIDV     residuals for this source,freq pair in Re, Im format
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,NREC), EST(MDIM,NREC) 
      REAL*8, INTENT(IN) :: AZMOD, FREQ 
      REAL*4, INTENT(IN) :: WGHTS(MDIM,NREC), DTMAX
      INTEGER*4, INTENT(IN) :: IPERM_REC(NDIM*NREC), MDIM, NOBS2, NDIM,NREC, &
                               IRESTP, IBPHASE
      REAL*4, INTENT(OUT) :: RESIDV(NOBS2) 
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      COMPLEX*8 HDIFFER, QN, QE, QZ, UEST, VEST, WEST, CZERO, &
                CRESIDU, CRESIDV, CRESIDW, CRESIDN, CRESIDE, NEST, EEST, &
                HFACT
      REAL*4 SPHASE, DPHI
      INTEGER*4 INDX, JNDX, IREC
      LOGICAL*4 LUSE(3) 
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      COMPLEX*8 CRESID
      REAL*4, PARAMETER :: TWOPI = 6.283185307179586
      REAL*4, PARAMETER :: PI = 3.14159265358979323846
!
!----------------------------------------------------------------------------------------!
!
!.... set indicees
      IERR = 0
      INDX = 0
      JNDX = NDIM*NREC !number of real observations for this source, freq pair 
      HFACT = HDIFFER(IBPHASE,FREQ) !option for half-differentiator
      RESIDV(1:NOBS2) = 0.0
!
!.... loop on receivers 
      DO 1 IREC=1,NREC
!
!....... get the estimates and observations, estimates as complex numbers
         QN = CZERO
         QE = CZERO
         QZ = CZERO
         UEST = CZERO
         VEST = CZERO
         WEST = CZERO
         LUSE(1:3) = .FALSE.
         IF (CABS(OBS(1,IREC)) > 0.0) THEN
            QN = OBS(1,IREC)
            LUSE(1) = .TRUE.
         ENDIF
         IF (CABS(OBS(2,IREC)) > 0.0) THEN
            QE = OBS(2,IREC)
            LUSE(2) = .TRUE.
         ENDIF
         IF (CABS(OBS(3,IREC)) > 0.0) THEN
            QZ = OBS(3,IREC)
            LUSE(3) = .TRUE.
         ENDIF
         UEST = EST(1,IREC)
         VEST = EST(2,IREC)
         WEST = EST(3,IREC)
!
!....... counterclockwise rotation from (U,V,W) to (N,E,Z)
         CALL CROTATE(+SNGL(AZMOD),UEST,VEST, NEST,EEST)
!
!....... calculate resiudals
         CRESIDU = CZERO
         CRESIDV = CZERO
         CRESIDW = CZERO
         CRESIDN = CZERO
         CRESIDE = CZERO
         IF (LUSE(1)) CRESIDN = CRESID(IRESTP,QN,NEST)*HFACT !u
         IF (LUSE(2)) CRESIDE = CRESID(IRESTP,QE,EEST)*HFACT !v
         IF (LUSE(3)) CRESIDW = CRESID(IRESTP,QZ,WEST)*HFACT !z
         CALL CROTATE(-SNGL(AZMOD),CRESIDN,CRESIDE, CRESIDU,CRESIDV)
!
!....... rotate back
         IF (LUSE(1) .and. LUSE(2)) THEN
            !clockwise rotation back into (U,V,W) frame
            CALL CROTATE(-SNGL(AZMOD),CRESIDN,CRESIDE, CRESIDU,CRESIDV)
            INDX =             IPERM_REC( (IREC - 1)*NDIM + 1) !INDX + 1
            JNDX = NDIM*NREC + IPERM_REC( (IREC - 1)*NDIM + 1) !JNDX + 1
            RESIDV(INDX) = REAL(CRESIDU)*WGHTS(1,IREC)
            RESIDV(JNDX) = IMAG(CRESIDU)*WGHTS(1,IREC)
            IF (DTMAX > 0.0) THEN
               DPHI = SPHASE(QN) - SPHASE(NEST)
               IF (DPHI <-PI) DPHI = DPHI + TWOPI 
               IF (DPHI >+PI) DPHI = DPHI - TWOPI 
               IF (ABS(DPHI/(TWOPI*SNGL(FREQ))) > DTMAX) THEN
                  RESIDV(INDX) = 0.0
                  RESIDV(JNDX) = 0.0
               ENDIF
            ENDIF
            INDX =             IPERM_REC( (IREC - 1)*NDIM + 2) !INDX + 1
            JNDX = NDIM*NREC + IPERM_REC( (IREC - 1)*NDIM + 2) !JNDX + 1
            RESIDV(INDX) = REAL(CRESIDV)*WGHTS(2,IREC)
            RESIDV(JNDX) = IMAG(CRESIDV)*WGHTS(2,IREC)
            IF (DTMAX > 0.0) THEN
               DPHI = SPHASE(QE) - SPHASE(EEST) 
               IF (DPHI <-PI) DPHI = DPHI + TWOPI 
               IF (DPHI >+PI) DPHI = DPHI - TWOPI 
               IF (ABS(DPHI/(TWOPI*SNGL(FREQ))) > DTMAX) THEN
                  RESIDV(INDX) = 0.0
                  RESIDV(JNDX) = 0.0
               ENDIF
            ENDIF
         ENDIF
!....... w
         INDX =             IPERM_REC( (IREC - 1)*NDIM + 3) !INDX + 1
         JNDX = NDIM*NREC + IPERM_REC( (IREC - 1)*NDIM + 3) !JNDX + 1
         IF (LUSE(3)) CRESIDW = CRESID(IRESTP,QZ,WEST)*HFACT
         RESIDV(INDX) = REAL(CRESIDW)*WGHTS(3,IREC)
         RESIDV(JNDX) = IMAG(CRESIDW)*WGHTS(3,IREC)
         IF (DTMAX > 0.0) THEN
            DPHI = SPHASE(QZ) - SPHASE(WEST)
            IF (DPHI <-PI) DPHI = DPHI + TWOPI 
            IF (DPHI >+PI) DPHI = DPHI - TWOPI 
            IF (ABS(DPHI/(TWOPI*SNGL(FREQ))) > DTMAX) THEN
               RESIDV(INDX) = 0.0
               RESIDV(JNDX) = 0.0
            ENDIF
         ENDIF
    1 CONTINUE 
!
!.... check indices 
!     IF (INDX /= NDIM*NREC) THEN
!        WRITE(*,*) 'oe2resv: Error indx /= ndim*nrec!',INDX,NDIM,NREC
!        IERR = 1
!     ENDIF 
!     IF (JNDX /= NOBS2) THEN
!        WRITE(*,*) 'oe2resv: Error jndx /= nobs!', JNDX, NOBS2
!        IERR = 2
!     ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SET_RCV_MOD(MDIM,MFREQ,MFREQ_LOC,              &
                             NFREQ, NDIM,NFREQ_REF,NREC,        &
                             LSURF, CFTYPE, RECLOC, RECV, IERR) 
!
!     Sets the updated receiver response function in reccloc onto the global 
!     RRF matrix.  This is for modeling frequencies
!
!     INPUT      MEANING
!     -----      ------- 
!     CFTYPE     frequency is for a body wave or surface wave
!     LSURF      True -> surface waves
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MFREQ_LOC  leading dimension
!     NFREQ      number of frequencies 
!     NFREQ_REF  number of inversion frequencies for surface or body waves
!     NREC       number of receivers 
!     RECLOC     RRFs for surface or body wave modeling
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     RECV       all RRFs 
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: RECLOC(MDIM,MFREQ_LOC,NREC) 
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MFREQ_LOC, NDIM, NFREQ, NFREQ_REF, NREC
      LOGICAL*4, INTENT(IN) :: LSURF
      COMPLEX*8, INTENT(INOUT) :: RECV(MDIM,MFREQ,NREC) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 IFREQ, JFREQ, IREC, I
      LOGICAL*4 LBODY, LCOPY
    
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0 
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies 
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         LCOPY = .FALSE.
         IF (LSURF) THEN
            IF (CFTYPE(IFREQ) == 'S') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ELSE
            IF (CFTYPE(IFREQ) == 'B') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ENDIF
         IF (LCOPY) THEN
            DO 2 IREC=1,NREC
               DO 3 I=1,NDIM
                  RECV(I,IFREQ,IREC) = RECLOC(I,JFREQ,IREC)
    3          CONTINUE !loop on components
    2       CONTINUE !loop on receivers 
         ENDIF !end check on copying
    1 CONTINUE !loop on frequencies
!
!.... check on sizes 
      IF (JFREQ /= NFREQ_REF) THEN
         WRITE(*,*) 'set_rcv_mod: Error jfreq /= nfreq_ref!'
         IERR = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_RCV_MOD(MDIM,MFREQ,MFREQ_LOC,              &
                             NFREQ, NDIM,NFREQ_REF,NREC,        &
                             LSURF, CFTYPE, RECV, RECLOC, IERR) 
!
!     Gets the RRF for modeling frequencies onto a local model
!
!     INPUT      MEANING
!     -----      ------- 
!     CFTYPE     frequency is for a body wave or surface wave
!     LSURF      True -> surface waves
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MFREQ_LOC  leading dimension
!     NFREQ      number of frequencies 
!     NFREQ_REF  number of inversion frequencies for surface or body waves
!     NREC       number of receivers
!     RECV       all RRFs
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     RECLOC     RRFs for surface or body wave modeling
!
!.... variable declarations
      IMPLICIT NONE 
      COMPLEX*8, INTENT(IN) :: RECV(MDIM,MFREQ,NREC) 
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MFREQ_LOC, NDIM, NFREQ, NFREQ_REF, NREC 
      LOGICAL*4, INTENT(IN) :: LSURF
      COMPLEX*8, INTENT(OUT) :: RECLOC(MDIM,MFREQ_LOC,NREC)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      INTEGER*4 IFREQ, JFREQ, IREC, I
      LOGICAL*4 LBODY, LCOPY
    
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies 
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         LCOPY = .FALSE.
         IF (LSURF) THEN
            IF (CFTYPE(IFREQ) == 'S') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ELSE
            IF (CFTYPE(IFREQ) == 'B') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ENDIF
         IF (LCOPY) THEN
            DO 2 IREC=1,NREC
               DO 3 I=1,NDIM
                  RECLOC(I,JFREQ,IREC) = RECV(I,IFREQ,IREC)
    3          CONTINUE !loop on components
    2       CONTINUE !loop on receivers 
         ENDIF !end check on copying
    1 CONTINUE !loop on frequencies
!
!.... check on sizes 
      IF (JFREQ /= NFREQ_REF) THEN
         WRITE(*,*) 'get_rcv_mod: Error jfreq /= nfreq_ref!'
         IERR = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SET_SRC_MOD(MFREQ,MFREQ_LOC, NFREQ,NSRC,  &
                             NFREQ_REF,NSRC_REF,  LSURF, CFTYPE, SRCTYP, &
                             SRCLOC, SOURCE, IERR) 
!
!     Sets the updated source time functions in srcloc onto the global 
!     STF matrix.  This is for modeling frequencies
!
!     INPUT      MEANING
!     -----      ------- 
!     CFTYPE     frequency is for a body wave or surface wave
!     LSURF      True -> surface waves
!     MFREQ      leading dimension
!     MFREQ_LOC  leading dimension
!     NFREQ      number of frequencies 
!     NFREQ_REF  number of inversion frequencies for surface or body waves
!     NSRC       number of sources
!     NSRC_REF   number of surface or body wave sources
!     SRCLOC     STFs for surface or body wave modeling 
!     SRCTYP     body or surface wave source
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     SOURCE     all STFs 
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: SRCLOC(MFREQ_LOC,NSRC_REF) 
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ, MFREQ_LOC, NFREQ, NSRC, NFREQ_REF, NSRC_REF 
      LOGICAL*4, INTENT(IN) :: LSURF
      COMPLEX*8, INTENT(INOUT) :: SOURCE(MFREQ,NSRC) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 IFREQ, JFREQ, ISRC, JSRC   
      LOGICAL*4 LBODY, LCOPY
    
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies 
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         LCOPY = .FALSE.
         IF (LSURF) THEN
            IF (CFTYPE(IFREQ) == 'S') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ELSE
            IF (CFTYPE(IFREQ) == 'B') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ENDIF
         IF (LCOPY) THEN
            JSRC = 0
            DO 2 ISRC=1,NSRC
               IF (LSURF .AND. SRCTYP(ISRC)(1:1) == 'S' .OR. &
                   LBODY .AND. SRCTYP(ISRC)(1:1) == 'P') THEN
                  JSRC = JSRC + 1
                  SOURCE(IFREQ,ISRC) = SRCLOC(JFREQ,JSRC)
               ENDIF
    2      CONTINUE
           IF (JSRC /= NSRC_REF) THEN
              WRITE(*,*) 'set_src_mod: Error jsrc /= nsrc_ref!'
              IERR = 1
              RETURN
           ENDIF
        ENDIF
    1 CONTINUE
!
!.... check on sizes 
      IF (JFREQ /= NFREQ_REF) THEN
         WRITE(*,*) 'set_src_mod: Error jfreq /= nfreq_ref!'
         IERR = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_SRC_MOD(MFREQ,MFREQ_LOC, NFREQ,NSRC,   &
                             NFREQ_REF,NSRC_REF,            &
                             LSURF, CFTYPE, SRCTYP, SOURCE, &
                             SRCLOC, IERR)
!
!     Gets the source time functions for modeling frequencies for surface or body waves
!
!     INPUT      MEANING
!     -----      ------- 
!     CFTYPE     frequency is for a body wave or surface wave
!     LSURF      True -> surface waves
!     MFREQ      leading dimension
!     MFREQ_LOC  leading dimension
!     NFREQ      number of frequencies 
!     NFREQ_REF  number of inversion frequencies for surface or body waves
!     NSRC       number of sources
!     NSRC_REF   number of surface or body wave sources
!     SOURCE     STFs for all frequencies sources 
!     SRCTYP     body or surface wave source
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     SRCLOC     local sources for only surface or body waves
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: SOURCE(MFREQ,NSRC)
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ, MFREQ_LOC, NFREQ, NSRC, NFREQ_REF, NSRC_REF
      LOGICAL*4, INTENT(IN) :: LSURF
      COMPLEX*8, INTENT(OUT) :: SRCLOC(MFREQ_LOC,NSRC_REF)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 IFREQ, JFREQ, ISRC, JSRC
      LOGICAL*4 LBODY, LCOPY

!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies 
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         LCOPY = .FALSE.
         IF (LSURF) THEN
            IF (CFTYPE(IFREQ) == 'S') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ELSE
            IF (CFTYPE(IFREQ) == 'B') THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ENDIF
         IF (LCOPY) THEN
            JSRC = 0
            DO 2 ISRC=1,NSRC
               IF (LSURF .AND. SRCTYP(ISRC)(1:1) == 'S' .OR. &
                   LBODY .AND. SRCTYP(ISRC)(1:1) == 'P') THEN
                  JSRC = JSRC + 1
                  SRCLOC(JFREQ,JSRC) = SOURCE(IFREQ,ISRC)
               ENDIF
    2      CONTINUE
           IF (JSRC /= NSRC_REF) THEN
              WRITE(*,*) 'get_src_mod: Error jsrc /= nsrc_ref!'
              IERR = 1
              RETURN
           ENDIF
        ENDIF
    1 CONTINUE
!
!.... check on sizes 
      IF (JFREQ /= NFREQ_REF) THEN
         WRITE(*,*) 'get_src_mod: Error jfreq /= nfreq_ref!'
         IERR = 1
      ENDIF
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SET_SRC_INV(MFREQ,MFREQ_LOC, NFREQ,NSRC,  &
                             NFREQ_REF,NSRC_REF,           &
                             LSURF, CFTYPE, SRCTYP, LINVF, SRCLOC, SOURCE, IERR) 
!
!     Sets the updated source time functions in srcloc onto the global 
!     STF matrix.  This is for inversion frequencies
!
!     INPUT      MEANING
!     -----      ------- 
!     CFTYPE     frequency is for a body wave or surface wave
!     LINVF      True -> inversion frequency
!     LSURF      True -> surface waves
!     MFREQ      leading dimension
!     MFREQ_LOC  leading dimension
!     NFREQ      number of frequencies 
!     NFREQ_REF  number of inversion frequencies for surface or body waves
!     NSRC       number of sources
!     NSRC_REF   number of surface or body wave sources
!     SRCLOC     local sources 
!     SRCTYP     body or surface wave source
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     SOURCE     all STFs
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: SRCLOC(MFREQ_LOC,*) 
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ, MFREQ_LOC, NFREQ, NSRC, NFREQ_REF, NSRC_REF 
      LOGICAL*4, INTENT(IN) :: LINVF(NFREQ), LSURF
      COMPLEX*8, INTENT(INOUT) :: SOURCE(MFREQ,*) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 IFREQ, JFREQ, ISRC, JSRC   
      LOGICAL*4 LBODY, LCOPY
      
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0 
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies 
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         LCOPY = .FALSE.
         IF (LSURF) THEN
            IF (CFTYPE(IFREQ) == 'S' .AND. LINVF(IFREQ)) THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ELSE
            IF (CFTYPE(IFREQ) == 'B' .AND. LINVF(IFREQ)) THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ENDIF
         IF (LCOPY) THEN
            JSRC = 0
            DO 2 ISRC=1,NSRC
               IF (LSURF .AND. SRCTYP(ISRC)(1:1) == 'S' .OR. &
                   LBODY .AND. SRCTYP(ISRC)(1:1) == 'P') THEN
                  JSRC = JSRC + 1
                  SOURCE(IFREQ,ISRC) = SRCLOC(JFREQ,JSRC)
               ENDIF
    2      CONTINUE
           IF (JSRC /= NSRC_REF) THEN
              WRITE(*,*) 'set_src_inv: Error jsrc /= nsrc_ref!'
              IERR = 1
              RETURN
           ENDIF
        ENDIF
    1 CONTINUE
!
!.... check on sizes 
      IF (JFREQ /= NFREQ_REF) THEN
         WRITE(*,*) 'set_src_inv: Error jfreq /= nfreq_ref!'
         IERR = 1 
      ENDIF
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_SRC_INV(MFREQ,MFREQ_LOC, NFREQ,NSRC,  &
                             NFREQ_REF,NSRC_REF, LSURF, &
                             CFTYPE, SRCTYP, LINVF, SOURCE, &
                             SRCLOC, IERR) 
!
!     Gets the source time functions for inversion frequencies  
!
!     INPUT      MEANING
!     -----      ------- 
!     CFTYPE     frequency is for a body wave or surface wave
!     LINVF      True -> inversion frequency
!     LSURF      True -> surface waves
!     MFREQ      leading dimension
!     MFREQ_LOC  leading dimension
!     NFREQ      number of frequencies 
!     NFREQ_REF  number of inversion frequencies for surface or body waves
!     NSRC       number of sources
!     NSRC_REF   number of surface or body wave sources
!     SOURCE     all STFs 
!     SRCTYP     body or surface wave source
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     SRCLOC     STFs for either surface or body waves inversion frequencies 
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: SOURCE(MFREQ,NSRC)  
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ, MFREQ_LOC, NFREQ, NSRC, NFREQ_REF, NSRC_REF 
      LOGICAL*4, INTENT(IN) :: LINVF(NFREQ), LSURF
      COMPLEX*8, INTENT(OUT) :: SRCLOC(MFREQ_LOC,NSRC_REF) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 IFREQ, JFREQ, ISRC, JSRC   
      LOGICAL*4 LBODY, LCOPY
    
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies 
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         LCOPY = .FALSE.
         IF (LSURF) THEN
            IF (CFTYPE(IFREQ) == 'S' .AND. LINVF(IFREQ)) THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ELSE
            IF (CFTYPE(IFREQ) == 'B' .AND. LINVF(IFREQ)) THEN
               JFREQ = JFREQ + 1
               LCOPY = .TRUE.
            ENDIF
         ENDIF
         IF (LCOPY) THEN
            JSRC = 0
            DO 2 ISRC=1,NSRC
               IF (LSURF .AND. SRCTYP(ISRC)(1:1) == 'S' .OR. &
                   LBODY .AND. SRCTYP(ISRC)(1:1) == 'P') THEN
                  JSRC = JSRC + 1
                  SRCLOC(JFREQ,JSRC) = SOURCE(IFREQ,ISRC) 
               ENDIF
    2      CONTINUE
           IF (JSRC /= NSRC_REF) THEN
              WRITE(*,*) 'get_src_inv: Error jsrc /= nsrc_ref!'
              IERR = 1
              RETURN
           ENDIF
        ENDIF
    1 CONTINUE
!
!.... check on sizes 
      IF (JFREQ /= NFREQ_REF) THEN
         WRITE(*,*) 'get_src_inv: Error jfreq /= nfreq_ref!'
         IERR = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_FREQ_INV(NFREQ, NFREQ_REF, LSURF, CFTYPE,LINVF,FREQ,  &
                              FREQ_LOC,IERR)
! 
!     Get the frequency list corresponding to the inversion frequencies
!
!     INPUT      MEANING
!     -----      -------
!     CFTYPE     S -> surface wave frequency, B -> body wave
!     FREQ       frequency list (Hz)
!     LINVF      True -> frequency is an inversion frequency
!     LSURF      True -> want surface wave frequencies
!     NFREQ      number of frequencies
!     NFREQ_REF  target number of frequencies, e.g. frq%nfreq_bdy_inv
!
!     OUTPUT     MEANING
!     ------     -------
!     FREQ_LOC   local frequency list for surface or body waves
!     IERR       error flag
!
!.... variable declarations  
      IMPLICIT NONE 
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      INTEGER*4, INTENT(IN) :: NFREQ, NFREQ_REF 
      LOGICAL*4, INTENT(IN) :: LINVF(NFREQ), LSURF
      REAL*8, INTENT(OUT) :: FREQ_LOC(NFREQ_REF)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      INTEGER*4 IFREQ, JFREQ 
      LOGICAL*4 LBODY 
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         IF (LSURF) THEN 
            IF (CFTYPE(IFREQ) == 'S' .AND. LINVF(IFREQ)) THEN 
               JFREQ = JFREQ + 1
               FREQ_LOC(JFREQ) = FREQ(IFREQ)
            ENDIF
         ELSE 
            IF (CFTYPE(IFREQ) == 'B' .AND. LINVF(IFREQ)) THEN 
               JFREQ = JFREQ + 1
               FREQ_LOC(JFREQ) = FREQ(IFREQ)
            ENDIF
         ENDIF
    1 CONTINUE
      IF (JFREQ /= NFREQ_REF) THEN 
         WRITE(*,*) 'get_freq_mod: Error jfreq /= nfreq_ref!'
         IERR = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_FREQ_MOD(NFREQ, NFREQ_REF, LSURF, CFTYPE,FREQ,  &
                              FREQ_LOC,IERR)
! 
!     Get the frequency list corresponding to the modeling frequencies
!
!     INPUT      MEANING
!     -----      -------
!     CFTYPE     S -> surface wave frequency, B -> body wave
!     FREQ       frequency list (Hz)
!     LSURF      True -> want surface wave frequencies
!     NFREQ      number of frequencies
!     NFREQ_REF  target number of frequencies, e.g. frq%nfreq_bdy
!
!     OUTPUT     MEANING
!     ------     -------
!     FREQ_LOC   local frequency list for surface or body waves
!     IERR       error flag
!
!.... variable declarations 
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      INTEGER*4, INTENT(IN) :: NFREQ, NFREQ_REF 
      LOGICAL*4, INTENT(IN) :: LSURF
      REAL*8, INTENT(OUT) :: FREQ_LOC(NFREQ_REF)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 IFREQ, JFREQ 
      LOGICAL*4 LBODY 
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface wave
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... loop on frequencies
      JFREQ = 0
      DO 1 IFREQ=1,NFREQ
         IF (LSURF) THEN
            IF (CFTYPE(IFREQ) == 'S') THEN
               JFREQ = JFREQ + 1
               FREQ_LOC(JFREQ) = FREQ(IFREQ)
            ENDIF
         ELSE
            IF (CFTYPE(IFREQ) == 'B') THEN
               JFREQ = JFREQ + 1
               FREQ_LOC(JFREQ) = FREQ(IFREQ)
            ENDIF
         ENDIF
    1 CONTINUE
      IF (JFREQ /= NFREQ_REF) THEN
         WRITE(*,*) 'get_freq_mod: Error jfreq /= nfreq_ref!'
         IERR = 1
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE COBJ_SB(LSURF,DT_MAX,AZMOD, INV,FRQ,SRC,RCV, OBJ, IERR) 
!
!     Calculates the objective function for surface or body waves
      IMPLICIT NONE 
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO)  INV  
      TYPE (SRC_INFO)  SRC  
      TYPE (FRQ_INFO)  FRQ  
      TYPE (RECV_INFO) RCV
      REAL*8, INTENT(IN) :: AZMOD
      REAL*4, INTENT(IN) :: DT_MAX
      LOGICAL*4, INTENT(IN) :: LSURF
      REAL*4, INTENT(OUT) :: OBJ
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      COMPLEX*8, ALLOCATABLE :: OBS_LOC(:,:,:,:), EST_LOC(:,:,:,:) 
      REAL*8, ALLOCATABLE :: FREQ_LOC(:)
      REAL*4, ALLOCATABLE :: WGHT_LOC(:,:,:,:)
      !COMPLEX*8 CPHM2CM 
      REAL*4 COBJ4_COVD, COBJ4!, PHASE, AMP 
      INTEGER*4 NFREQ, NREC, NSRC !,IFREQ,IREC,ISRC,I 
!     integer*4 ifreq,irec,isrc
!
!----------------------------------------------------------------------------------------!
!
!.... quick return 
      IERR = 0
      NREC = rcv%NREC
      IF (LSURF) THEN 
         IF (.NOT.inv%LSURF) THEN 
            WRITE(*,*) 'cobj_sb: We are not inverting surface waves, returning...'
            RETURN
         ENDIF
         NFREQ = frq%NFREQ_SRF_INV
         NSRC  = src%NSRC_SRF
      ELSE 
         IF (.NOT.inv%LBODY) THEN 
            WRITE(*,*) 'cobj_sb: We are not inverting body waves, returning...'
            RETURN
         ENDIF
         NFREQ = frq%NFREQ_BDY_INV
         NSRC  = src%NSRC_BDY
      ENDIF
!
!.... set space and extract 
      ALLOCATE( OBS_LOC(NDIM,NFREQ,rcv%NREC,NSRC))
      ALLOCATE( EST_LOC(NDIM,NFREQ,rcv%NREC,NSRC))
      ALLOCATE(WGHT_LOC(NDIM,NFREQ,rcv%NREC,NSRC))
      ALLOCATE(FREQ_LOC(NFREQ)) 
      CALL GET_FREQ_INV(frq%NFREQ, NFREQ, LSURF, frq%CFTYPE,frq%LINVF,frq%FREQ, &
                        FREQ_LOC,IERR)
      IF (LSURF) THEN
         CALL GET_OEINV_SB( .TRUE., NDIM,NFREQ,NREC, INV,FRQ,SRC,RCV, &
                           WGHT_LOC,EST_LOC,OBS_LOC,IERR)
      ELSE
         CALL GET_OEINV_SB(.FALSE., NDIM,NFREQ,NREC, INV,FRQ,SRC,RCV, &
                           WGHT_LOC,EST_LOC,OBS_LOC,IERR)
      ENDIF
      IF (IERR /= 0) THEN
         WRITE(*,*) 'cobj_sb: Error calling get_oeinv_sb'
         RETURN
      ENDIF
!
!.... if the data is unwrapped then we need to fix that
!     IF (inv%LUNWRAP) THEN
!        DO 1 IFREQ=1,NFREQ
!           DO 2 ISRC=1,NSRC 
!              DO 3 IREC=1,NREC
!                 DO 4 I=1,NDIM
!                    AMP   = REAL(EST_LOC(I,IFREQ,IREC,ISRC))
!                    PHASE = IMAG(EST_LOC(I,IFREQ,IREC,ISRC))
!                    EST_LOC(I,IFREQ,IREC,ISRC) = CPHM2CM(AMP,PHASE) 
!                    AMP   = REAL(OBS_LOC(I,IFREQ,IREC,ISRC))
!                    PHASE = IMAG(OBS_LOC(I,IFREQ,IREC,ISRC))
!                    OBS_LOC(I,IFREQ,IREC,ISRC) = CPHM2CM(AMP,PHASE)
!   4             CONTINUE
!   3          CONTINUE    
!   2       CONTINUE
!   1    CONTINUE
!     ENDIF
!
!.... rotate the estimates (u,v,w) -> (N,E,Z)
!     CALL ROTEST(NDIM,NFREQ,NREC, NDIM,NFREQ,NREC,NSRC,  &
!                 AZMOD,EST_LOC)
!     do ifreq=1,nfreq
!        do isrc=1,nsrc
!           do irec=1,nrec
!              if (lsurf) then
!              write(55,*) ifreq,irec,isrc, obs_loc(:,ifreq,irec,isrc)
!              write(55,*) ifreq,irec,isrc, est_loc(:,ifreq,irec,isrc)
!              else
!              write(56,*) ifreq,irec,isrc, obs_loc(:,ifreq,irec,isrc)
!              write(56,*) ifreq,irec,isrc, est_loc(:,ifreq,irec,isrc)
!              endif
!           enddo
!        enddo
!     enddo
!
!.... calculate objective function
      IF (     LSURF .AND. inv%LCOVD_SRF .OR. .NOT.LSURF .AND. inv%LCOVD_BDY) THEN
         OBJ = COBJ4_COVD(NDIM,NFREQ,NREC, NDIM,NFREQ,NREC,NSRC,                   &
                          inv%IRESTP,inv%IBPHASE, .TRUE.,inv%LDWGHT,               &
                          DT_MAX,AZMOD,inv%DX, FREQ_LOC, WGHT_LOC,OBS_LOC,EST_LOC, &
                          IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'cobj_sb: Error calling cobj4_covd!'
            RETURN
         ENDIF
      ELSE
!
!....... rotate the estimates (u,v,w) -> (N,E,Z)
         CALL ROTEST(NDIM,NFREQ,NREC, NDIM,NFREQ,NREC,NSRC,  &
                     AZMOD,EST_LOC)
         OBJ = COBJ4(NDIM,NFREQ,NREC, NFREQ,NREC,NSRC, NDIM,     &
                     .FALSE.,2,inv%IRESTP, 0.2,0.2,DT_MAX,       &
                     FREQ_LOC,WGHT_LOC,OBS_LOC,EST_LOC)
      ENDIF
!
!.... free space
      DEALLOCATE(OBS_LOC)
      DEALLOCATE(EST_LOC)
      DEALLOCATE(WGHT_LOC)
      DEALLOCATE(FREQ_LOC)
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE MASK_LIFREQ(NFREQ,LSURF,JOB, CFTYPE, LINVF,LFSAVE)
!
!     To generate search directions for the surface and body waves we calculate 
!     independent gradients/search directions that satisfy a line search.  
!     So we aren't re-doubling our efforts on generating information we won't
!     use we mask the surface wave inversion frequencies for a body wave 
!     inversion section and vice versa.
!
!     INPUT      MEANING
!     -----      ------- 
!     JOB        = 1 apply mask; otherwise restore linvf 
!     LSURF      True -> keeping surface waves; False -> keeping body waves 
!     NFREQ      number of total frequencies
!  
!     OUTPUT     MEANING
!     ------     -------
!     LFSAVE     job = 1 -> saved original copy of linvf
!                job = 2 -> used to restore original linvf
!     LINVF      job = 1 -> inversion frequencies for surface or body waves only
!                job = 2 -> restored to original
!
!.... variable declarations
      LOGICAL*4, INTENT(INOUT) :: LINVF(NFREQ), LFSAVE(NFREQ)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: NFREQ, JOB
      LOGICAL*4, INTENT(IN) :: LSURF
      INTEGER*4 IFREQ
!
!----------------------------------------------------------------------------------------!
!
!.... masking frequencies
      IF (JOB == 1) THEN
         DO 1 IFREQ=1,NFREQ
            LFSAVE(IFREQ) = LINVF(IFREQ)
            IF (LSURF) THEN !surface wave inversion section 
               IF (CFTYPE(IFREQ) == 'B') LINVF(IFREQ) = .FALSE. 
            ELSE !body wave inversion section
               IF (CFTYPE(IFREQ) == 'S') LINVF(IFREQ) = .FALSE.
            ENDIF
    1    CONTINUE
      ELSE !restoring frequencies
         LINVF(1:NFREQ) = LFSAVE(1:NFREQ)
      ENDIF
      RETURN
      END
