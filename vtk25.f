!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_ELMOD_VTK(PROJNM,MGNOD,MNPG,NNPG, NELEM, LISISO, 
     ;                          IENG,XLOCS,ZLOCS, DENS,ECOEFF)
! 
!     Plots the elastic model.  The output is a binary .vtk file 
! 
!     INPUT      MEANING
!     -----      ------- 
!     DENS       density at anchor nodes
!     ECOEFF     elastic coefficients at anchor nodes
!     IENG       connectivity vector for anchor nodes
!     LISISO     True -> elements are isotropic, i.e. lambda/mu
!                False -> elements are anisotropic
!     MGNOD      max number of anchor nodes
!     NELEM      number of elements 
!     NNPG       number of anchor nodes 
!     PROJNM     project name
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!    
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) ::  ECOEFF(MNPG,*), DENS(NNPG), 
     ;                       XLOCS(NNPG), ZLOCS(NNPG)
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MGNOD,MNPG,NNPG,NELEM 
      LOGICAL*4, INTENT(IN) :: LISISO 
!.... local variables
      CHARACTER(80) FILENM
      CHARACTER(1), ALLOCATABLE :: FPASS(:)
      REAL*8 ALPHA2, BETA2 
      REAL*4, ALLOCATABLE :: XPTS(:), ZPTS(:), DAT(:) 
      INTEGER*4, ALLOCATABLE :: IENGV(:), IENGNOD(:)
      LOGICAL*4 LEX, LISDIR 
! 
!----------------------------------------------------------------------!
! 
!.... extract mesh partition and convert IENG to a vector 
      NWORK = 4*NELEM          !worst case is all linear quads 
      ALLOCATE(IENGV(NWORK))   !will hold vectorized IENG vector 
      ALLOCATE(IENGNOD(NELEM)) !eventually will be input
      ALLOCATE(XPTS(NNPG))     !VisIt uses single precision
      ALLOCATE(ZPTS(NNPG))     !VisIt uses single precision
      NIENGV = 0
      DO 1 IELEM=1,NELEM
         IENGNOD(IELEM) = 4 !default to quads
         DO 2 IA=1,4 !NGNOD  
            INPG = IENG(IA,IELEM)
            IF (INPG.EQ.0) THEN
               IENGNOD(IELEM) = 3
               GOTO 20
            ENDIF
            NIENGV = NIENGV + 1
            IENGV(NIENGV) = INPG
    2    CONTINUE
   20    CONTINUE
    1 CONTINUE
! 
!.... copy over nodal locations
      IF (LISISO) THEN
         ALLOCATE(DAT(3*NNPG))  
      ELSE
         WRITE(*,*) 'plot_model_vtk: Error program lisiso = false'
      ENDIF
      DO 3 INPG=1,NNPG
         XPTS(INPG) = SNGL(XLOCS(INPG))
         ZPTS(INPG) = SNGL(ZLOCS(INPG))
         IF (LISISO) THEN
            LOC1 = INPG
            LOC2 = NNPG + INPG
            LOC3 = 2*NNPG + INPG 
            ALPHA2 = (ECOEFF(INPG,1) + 2.D0*ECOEFF(INPG,2))/DENS(INPG)
            BETA2 = ECOEFF(INPG,2)/DENS(INPG) 
            DAT(LOC1) = SNGL(DSQRT(ALPHA2))
            DAT(LOC2) = SNGL(DSQRT(BETA2))  
            DAT(LOC3) = SNGL(DENS(INPG)) 
         ELSE
            WRITE(*,*) 'plot_model_vtk: Error no anistropic ordering'
            RETURN
         ENDIF 
    3 CONTINUE
! 
!.... write vtk file
      LEX = LISDIR('./model_vtk')
!#ifdef INTEL
!      INQUIRE(DIRECTORY='./model_vtk',EXIST=LEX)
!#else
!      INQUIRE(FILE='./model_vtk',EXIST=LEX)
!#endif
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./model_vtk')
      CALL NULLS(80,FILENM)
      FILENM = './model_vtk/'//TRIM(PROJNM)//'_model.vtk'
      FILENM = ADJUSTL(FILENM) 
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL)) 
      DO 4 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    4 CONTINUE
      CALL VTK_MODEL(FPASS,LENFL, NNPG,NNPG,NELEM,NIENGV,
     ;               IENGNOD,IENGV, XPTS,ZPTS,DAT)
      DEALLOCATE(IENGV)
      DEALLOCATE(IENGNOD)
      DEALLOCATE(XPTS)
      DEALLOCATE(ZPTS)
      DEALLOCATE(DAT) 
      DEALLOCATE(FPASS) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DUMP_CSRFULL(PROJNM, N,NZERO,PART,IRPTR,JCPTR) 
!
!     Writes the full compressed sparse row storage for plotting
!
!     INPUT      MEANING
!     -----      -------
!     IRPTR      CSR row pointer
!     JCPTR      CSR column pointer
!     N          number of rows
!     NZERO      number of non-zeros
!     PART       DOFs partition number 
!     PROJNM     project name
!
!.... local variables
      CHARACTER(*), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: PART(N), IRPTR(N+1), JCPTR(NZERO), 
     ;                         N, NZERO 
      CHARACTER(80) FILENM
      INTEGER*4 I,J, J1,J2
      LOGICAL*4 LISDIR, LEX
      INTEGER*4, PARAMETER :: IUNIT = 55
!----------------------------------------------------------------------!
      LEX = LISDIR('./mesh_partition') 
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./mesh_partition')
      FILENM(1:80) = ' '
      FILENM = './mesh_partition/'//TRIM(ADJUSTL(PROJNM))//'_CSR.dat'
      FILENM = TRIM(ADJUSTL(FILENM))
      OPEN(UNIT=IUNIT,FILE=TRIM(ADJUSTL(FILENM)))
      WRITE(IUNIT,*) N,NZERO
      WRITE(IUNIT,*) (IRPTR(I),I=1,N+1)
      WRITE(IUNIT,*) (PART(I),I=1,N)
      DO 1 I=1,N
         J1 = IRPTR(I)
         J2 = IRPTR(I+1) - 1
         WRITE(IUNIT,*) (JCPTR(J),J=J1,J2)
    1 CONTINUE
      CLOSE(IUNIT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_ELCRESP_VTK(MGNOD,MDIM,MEN, NDIM,NNPG,NDOF,
     ;                            NLXI,NLETA,PROJNM, NELEM,NCOMP, ISRC, 
     ;                            IFTYPE,FREQ, LM,IENG, XLOCS,ZLOCS, U) 

      CHARACTER(*), INTENT(IN) :: PROJNM 
      COMPLEX*8, INTENT(IN) :: U(NDOF) 
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), FREQ  
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*)
!.... local variables
      CHARACTER(1), ALLOCATABLE :: FPASS(:) 
      REAL*4, ALLOCATABLE :: RMAG(:), PHASE(:), XPTS(:), ZPTS(:) 
      INTEGER*4, ALLOCATABLE :: IDPTR(:,:), IENGNOD(:), IENGV(:) 
      CHARACTER(80) FILENM 
      CHARACTER(12) CFREQ
      CHARACTER(10) PRE
      CHARACTER(5) CSRC, APP 
      REAL*4 PI180 
      LOGICAL*4 LEX, LISDIR
      PARAMETER(PI180 = 0.017453292519943295) 
      PARAMETER(LNULL =-5)
! 
!----------------------------------------------------------------------!
! 
!.... extract mesh partition and convert IENG to a vector
      NWORK = 4*NELEM            !worst case is all linear quads 
      ALLOCATE(IENGV(NWORK))     !will hold vectorized IENG vector 
      ALLOCATE(IENGNOD(NELEM))   !eventually will be input 
      LDD = NNPG 
      NIENGV = 0 
      DO 1 IELEM=1,NELEM
         IENGNOD(IELEM) = 4 !default to quad
         DO 2 IA=1,4 !NGNOD
            INPG = IENG(IA,IELEM)
            IF (INPG.EQ.0) THEN
               IENGNOD(IELEM) = 3 
               GOTO 20
            ENDIF
            NIENGV = NIENGV + 1 
            IENGV(NIENGV) = INPG 
    2    CONTINUE
   20    CONTINUE 
    1 CONTINUE
!
!.... set space 
      ALLOCATE(XPTS(NNPG))
      ALLOCATE(ZPTS(NNPG))
      ALLOCATE(RMAG (NCOMP*LDD)) !holds u and w magnitudes 
      ALLOCATE(PHASE(NCOMP*LDD)) !holds u and w phases
      ALLOCATE(IDPTR(NNPG,NDIM))
      CALL GENIDPTR(NNPG,MGNOD,MEN,MDIM, NNPG,NELEM,NDIM, 
     ;              NLXI,NLETA,LNULL, IENGNOD,IENG,LM, IDPTR) 
      DO 3 INPG=1,NNPG
         XPTS(INPG) = SNGL(XLOCS(INPG))
         ZPTS(INPG) = SNGL(ZLOCS(INPG))
         DO 4 I=1,NDIM
            INDX = (I - 1)*LDD + INPG
            IDOF = IDPTR(INPG,I)
            IF (IDOF.EQ.LNULL) THEN
               WRITE(*,*) 'plot_eqwave_vtk: Error calling genidptr'
               IERR = 1 
            ENDIF
            IF (IDOF.EQ.0) THEN
               RMAG(INDX) = 0.0 
               PHASE(INDX) = 0.0 
            ELSE
               RMAG(INDX) = CABS(U(IDOF))
               IF (RMAG(INDX).LT.1.11E-7) THEN
                  PHASE(INDX) = 0.0
               ELSE
                  PHASE(INDX) = ATAN2(IMAG(U(IDOF)),REAL(U(IDOF)))/PI180
               ENDIF
            ENDIF
    4    CONTINUE
    3 CONTINUE
      DEALLOCATE(IDPTR)

      CALL NULLS(80,FILENM)
      CALL NULLS(12,CFREQ)
      CALL NULLS(5,CSRC) 
      WRITE(CFREQ,'(F12.5)') FREQ
      CFREQ = ADJUSTL(CFREQ) 
      WRITE(CSRC,'(I5)') ISRC 
      CSRC = ADJUSTL(CSRC) 
      IF (IFTYPE.EQ.1) THEN
         PRE = './blk_vtk/'
         APP = '_blk-'
      ELSEIF (IFTYPE.EQ.2) THEN
         PRE = './frc_vtk/'
         APP = '_frc-'
      ELSE
         PRE = './wav_vtk/'
         APP = '_wav-'
      ENDIF
!
!.... this should not be done here b/c MPI, this is just a backup
      LEX = LISDIR(TRIM(PRE))
!#ifdef INTEL
!      INQUIRE(DIRECTORY=TRIM(PRE),EXIST=LEX)
!#else
!      INQUIRE(FILE=TRIM(PRE),EXIST=LEX)
!#endif
      IF (.NOT.LEX) THEN
         IF (IFTYPE.EQ.1) THEN
            CALL SYSTEM('mkdir ./blk_vtk')
         ELSEIF (IFTYPE.EQ.2) THEN
            CALL SYSTEM('mkdir ./frc_vtk')
         ELSE
            CALL SYSTEM('mkdir ./wav_vtk')
         ENDIF
      ENDIF
      FILENM = PRE//TRIM(ADJUSTL(PROJNM))//APP//TRIM(CFREQ)//'-'//
     ;         TRIM(CSRC)//'.vtk'
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL)) 
      DO 5 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    5 CONTINUE 
      CALL VTK_QRESP25D_EL(FPASS,LENFL, LDD,NNPG,NELEM,NIENGV, 
     ;                     IENGNOD,IENGV, XPTS,ZPTS,RMAG,PHASE) 
      DEALLOCATE(XPTS)
      DEALLOCATE(ZPTS) 
      DEALLOCATE(RMAG)
      DEALLOCATE(PHASE)
      DEALLOCATE(IENGV) 
      DEALLOCATE(IENGNOD) 
      DEALLOCATE(FPASS) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_BPRESID_VTK(MGNOD,MDIM,MEN, NDIM,NNPG,NDOF,
     ;                            NLXI,NLETA,PROJNM, NELEM,NCOMP, ISRC, 
     ;                            IREC,FREQ, LM,IENG, XLOCS,ZLOCS, U) 

      CHARACTER(*), INTENT(IN) :: PROJNM 
      COMPLEX*8, INTENT(IN) :: U(NDOF) 
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), FREQ  
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*)
!.... local variables
      CHARACTER(1), ALLOCATABLE :: FPASS(:) 
      REAL*4, ALLOCATABLE :: RMAG(:), PHASE(:), XPTS(:), ZPTS(:) 
      INTEGER*4, ALLOCATABLE :: IDPTR(:,:), IENGNOD(:), IENGV(:) 
      CHARACTER(80) FILENM 
      CHARACTER(12) CFREQ
      CHARACTER(5) CSRC, CREC 
      REAL*4 PI180 
      LOGICAL*4 LEX, LISDIR
      PARAMETER(PI180 = 0.017453292519943295) 
      PARAMETER(LNULL =-5)
! 
!----------------------------------------------------------------------!
! 
!.... extract mesh partition and convert IENG to a vector
      NWORK = 4*NELEM            !worst case is all linear quads 
      ALLOCATE(IENGV(NWORK))     !will hold vectorized IENG vector 
      ALLOCATE(IENGNOD(NELEM))   !eventually will be input 
      LDD = NNPG 
      NIENGV = 0 
      DO 1 IELEM=1,NELEM
         IENGNOD(IELEM) = 4 !default to quad
         DO 2 IA=1,4 !NGNOD
            INPG = IENG(IA,IELEM)
            IF (INPG.EQ.0) THEN
               IENGNOD(IELEM) = 3 
               GOTO 20
            ENDIF
            NIENGV = NIENGV + 1 
            IENGV(NIENGV) = INPG 
    2    CONTINUE
   20    CONTINUE 
    1 CONTINUE
!
!.... set space 
      ALLOCATE(XPTS(NNPG))
      ALLOCATE(ZPTS(NNPG))
      ALLOCATE(RMAG (NCOMP*LDD)) !holds u and w magnitudes 
      ALLOCATE(PHASE(NCOMP*LDD)) !holds u and w phases
      ALLOCATE(IDPTR(NNPG,NDIM))
      CALL GENIDPTR(NNPG,MGNOD,MEN,MDIM, NNPG,NELEM,NDIM,
     ;              NLXI,NLETA,LNULL, IENGNOD,IENG,LM, IDPTR)
      DO 3 INPG=1,NNPG
         XPTS(INPG) = SNGL(XLOCS(INPG))
         ZPTS(INPG) = SNGL(ZLOCS(INPG))
         DO 4 I=1,NDIM
            INDX = (I - 1)*LDD + INPG
            IDOF = IDPTR(INPG,I)
            IF (IDOF.EQ.LNULL) THEN
               WRITE(*,*) 'plot_eqwave_vtk: Error calling genidptr'
               IERR = 1
            ENDIF
            IF (IDOF.EQ.0) THEN
               RMAG(INDX) = 0.0
               PHASE(INDX) = 0.0
            ELSE
               RMAG(INDX) = CABS(U(IDOF))
               IF (RMAG(INDX).LT.1.11E-7) THEN
                  PHASE(INDX) = 0.0
               ELSE
                  PHASE(INDX) = ATAN2(IMAG(U(IDOF)),REAL(U(IDOF)))/PI180
               ENDIF
            ENDIF
    4    CONTINUE
    3 CONTINUE
      DEALLOCATE(IDPTR)
      CALL NULLS(80,FILENM)
      CALL NULLS(12,CFREQ)
      CALL NULLS(5,CSRC)
      CALL NULLS(5,CREC) 
      WRITE(CFREQ,'(F12.5)') FREQ
      CFREQ = ADJUSTL(CFREQ)
      WRITE(CSRC,'(I5)') ISRC
      CSRC = ADJUSTL(CSRC)
      WRITE(CREC,'(I5)') IREC 
      CREC = ADJUSTL(CREC) 
      LEX = LISDIR('./baw_vtk')
!#ifdef INTEL
!      INQUIRE(DIRECTORY=TRIM('./baw_vtk'),EXIST=LEX)
!#else
!      INQUIRE(FILE=TRIM('./baw_vtk'),EXIST=LEX)
!#endif
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./baw_vtk')
      FILENM = './baw_vtk/'//TRIM(ADJUSTL(PROJNM))//'_baw-'//
     ;         TRIM(CFREQ)//'-'//TRIM(CSRC)//'-'//TRIM(CREC)//'.vtk'
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 5 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    5 CONTINUE
      CALL VTK_QRESP25D_EL(FPASS,LENFL, LDD,NNPG,NELEM,NIENGV,
     ;                     IENGNOD,IENGV, XPTS,ZPTS,RMAG,PHASE)
      DEALLOCATE(XPTS)
      DEALLOCATE(ZPTS)
      DEALLOCATE(RMAG)
      DEALLOCATE(PHASE)
      DEALLOCATE(IENGV)
      DEALLOCATE(IENGNOD)
      DEALLOCATE(FPASS)
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_ELRESP_VTK(MGNOD,MNPG,NNPG,PROJNM, NELEM,ITER, 
     ;                           ISRC,NCOMP,IENG, XLOCS,ZLOCS,RESP) 
! 
!     INPUT      MEANING
!     -----      ------- 
!     IENG       anchor node pointer
!     ISRC       source number
!     ITER       time slice number
!     MGNOD      max number of anchor nodes
!     MNPG       max number of anchor node points
!     NCOMP      number of components in response, probably 2
!     NELEM      number of elements
!     NNPG       number of anchor nodes points
!     PROJNM     project name
!     RESP       time domain displacement 
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes 
!
!.... variable declarations 
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*4 XLOCS(NNPG), ZLOCS(NNPG), RESP(MNPG,NCOMP)
      INTEGER*4 IENG(MGNOD,*) 
!.... local variables 
      CHARACTER(80) FILENM
      CHARACTER(1), ALLOCATABLE :: FPASS(:)
      CHARACTER(5) CITER,CSRC
      REAL*4, ALLOCATABLE :: RESPV(:) 
      INTEGER*4, ALLOCATABLE :: IENGV(:),IENGNOD(:) 
      LOGICAL*4 LEX, LISDIR
! 
!----------------------------------------------------------------------!
! 
!.... extract mesh partition and convert IENG to a vector
      NWORK = 4*NELEM            !worst case is all linear quads 
      ALLOCATE(IENGV(NWORK))     !will hold vectorized IENG vector 
      ALLOCATE(IENGNOD(NELEM))   !eventually will be input 
      LDD = NNPG 
      ALLOCATE(RESPV(NCOMP*LDD)) !holds u and w responses 
      NIENGV = 0
      DO 1 IELEM=1,NELEM
         IENGNOD(IELEM) = 4 !default to quad
         DO 2 IA=1,4 !NGNOD
            INPG = IENG(IA,IELEM)
            IF (INPG.EQ.0) THEN
               IENGNOD(IELEM) = 3
               GOTO 20
            ENDIF
            NIENGV = NIENGV + 1
            IENGV(NIENGV) = INPG 
    2    CONTINUE
   20    CONTINUE 
    1 CONTINUE
! 
!.... turn response into a vector
      RESPV(1:NCOMP*LDD) = 0.0
      DO 3 INPG=1,NNPG
         DO 4 I=1,NCOMP
            ILOC = (I - 1)*LDD + INPG 
            RESPV(ILOC) = RESP(INPG,I)
    4    CONTINUE
    3 CONTINUE 
! 
!.... check on directory (shouldn't be here b/c MPI) 
      LEX = LISDIR('./movie_vtk')
!#ifdef INTEL
!      INQUIRE(DIRECTORY='./movie_vtk',EXIST=LEX)
!#else
!      INQUIRE(FILE='./movie_vtk',EXIST=LEX)
!#endif
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./movie_vtk')
!
!.... write vtk file
      CALL NULLS(5,CSRC)
      WRITE(CSRC,'(I5)') ISRC
      CALL NULLS(5,CITER)
      WRITE(CITER,'(I5)') ITER
      CALL NULLS(80,FILENM)
      FILENM = './movie_vtk/'//TRIM(PROJNM)//'-'//
     ;         TRIM(ADJUSTL(CSRC) )//'-'//
     ;         TRIM(ADJUSTL(CITER))//'_ts.vtk'
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 5 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    5 CONTINUE
      CALL VTK_RRESP25D_EL(FPASS,LENFL, LDD,NNPG,NELEM,NIENGV,
     ;                     IENGNOD,IENGV, RESPV,XLOCS,ZLOCS)
      DEALLOCATE(IENGV)
      DEALLOCATE(IENGNOD)
      DEALLOCATE(RESPV)
      DEALLOCATE(FPASS) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_MESHP_VTK(PROJNM,MGNOD,MDIM,MEN, NNPG,NELEM,NDOF,
     ;                          NDIM,NEN, IDPART,IENG,LM, XLOCS,ZLOCS)
! 
!     Plots the mesh partition.  The output is a binary .vtk file 
! 
!     INPUT      MEANING
!     -----      ------- 
!     IENG       connectivity vector for anchor nodes
!     IDPART     holds partition number of each DOF 
!     ISABS      determines of element is absorbing or not
!     MGNOD      max number of anchor nodes
!     NABS       number of absorbing elements
!     NELEM      number of elements 
!     NNPG       number of anchor nodes 
!     PROJNM     project name
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!    
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: XLOCS(NNPG),ZLOCS(NNPG)
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), 
     ;           IDPART(NDOF), MGNOD,MDIM,MEN, NNPG,NELEM,NDOF,NDIM,NEN
!.... local variables
      CHARACTER(80) FILENM
      CHARACTER(1), ALLOCATABLE :: FPASS(:) 
      REAL*4, ALLOCATABLE :: XPTS(:), ZPTS(:), PART(:)
      INTEGER*4, ALLOCATABLE :: IENGV(:), IENGNOD(:), MYHIT(:)
      LOGICAL*4 LEX, LISDIR
! 
!----------------------------------------------------------------------!
! 
!.... extract mesh partition and convert IENG to a vector 
      NWORK = 4*NELEM          !worst case is all linear quads 
      ALLOCATE(IENGV(NWORK))   !will hold vectorized IENG vector 
      ALLOCATE(PART(NELEM))    !will hold entire element partition
      ALLOCATE(IENGNOD(NELEM)) !eventually will be input
      ALLOCATE(XPTS(NNPG))     !VisIt uses single precision
      ALLOCATE(ZPTS(NNPG))     !VisIt uses single precision
      I1 = 0
      J1 = 0
      NIENGV = 0
      NPARTS = MAXVAL(IDPART)
      ALLOCATE(MYHIT(NPARTS))
      DO 1 IELEM=1,NELEM
         NAVG = 0
         MYHIT(1:NPARTS) = 0
         DO 2 IAE=1,NEN
            DO 3 I=1,NDIM 
               IDOF = LM(I,IAE,IELEM)
               IF (IDOF.GT.0) THEN !dof
                  INDX = IDPART(IDOF)
                  MYHIT(INDX) = MYHIT(INDX) + 1
               ENDIF
    3       CONTINUE
    2    CONTINUE 
         IRES = MAXLOC(MYHIT,1) 
         IF (IRES.GT.0) THEN
            PART(IELEM) = FLOAT(IRES) 
         ELSE
            PART(IELEM) = 0.0
         ENDIF
         IENGNOD(IELEM) = 4 !default to quads
         DO 4 IA=1,4 !NGNOD  
            INPG = IENG(IA,IELEM)
            IF (INPG.EQ.0) THEN
               IENGNOD(IELEM) = 3 
               GOTO 40
            ENDIF 
            NIENGV = NIENGV + 1
            IENGV(NIENGV) = INPG
    4    CONTINUE
   40    CONTINUE 
    1 CONTINUE
      DEALLOCATE(MYHIT)
! 
!.... copy over nodal locations
      DO 5 INPG=1,NNPG
         XPTS(INPG) = SNGL(XLOCS(INPG))
         ZPTS(INPG) = SNGL(ZLOCS(INPG))
    5 CONTINUE
! 
!.... write vtk file
      LEX = LISDIR('./mesh_partition')
!#ifdef INTEL
!      INQUIRE(DIRECTORY='./mesh_partition',EXIST=LEX)
!#else
!      INQUIRE(FILE='./mesh_partition',EXIST=LEX)
!#endif
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./mesh_partition')
      CALL NULLS(80,FILENM)
      CALL DBLANK(80,PROJNM,LENFL)
      FILENM = './mesh_partition/'//PROJNM(1:LENFL)//'_mesh_part.vtk'
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 6 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    6 CONTINUE
      CALL VTK_MESHP2D(FPASS,LENFL, NNPG,NELEM,NIENGV,
     ;                 IENGNOD,IENGV,PART, XPTS,ZPTS)
      DEALLOCATE(IENGV)
      DEALLOCATE(PART)
      DEALLOCATE(IENGNOD)
      DEALLOCATE(XPTS)
      DEALLOCATE(ZPTS)
      DEALLOCATE(FPASS) 
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GENIDPTR(MNPG,MGNOD,MEN,MDIM, NNPG,NELEM,NDIM,
     ;                    NLXI,NLETA,LNULL, IENGNOD,IENG,LM, IDPTR)
      INTEGER*4 LM(MDIM,MEN,NELEM), IENGNOD(NELEM), IENG(MGNOD,NELEM)
      INTEGER*4 IDPTR(MNPG,NDIM)
! 
!----------------------------------------------------------------------!
!
      IDPTR(1:NNPG,1:NDIM) = LNULL
      DO 1 IELEM=1,NELEM
         IF (IENGNOD(IELEM).EQ.4) THEN
            NGNOD = 4
         ELSE
            NGNOD = 3
         ENDIF
         DO 2 IA=1,NGNOD
            ILOC = IENG(IA,IELEM)
            IF (ILOC.GT.0) THEN
               DO 3 I=1,NDIM
                  IF (IDPTR(ILOC,I).EQ.LNULL) THEN
                     IAE = IA2IAE(NLXI,NLETA,IA)
                     IDOF = LM(I,IAE,IELEM)
                     IF (IDOF.EQ.0) THEN
                        IDPTR(ILOC,I) = 0
                     ELSE
                        IDPTR(ILOC,I) = IDOF
                     ENDIF
                  ENDIF
    3          CONTINUE
            ENDIF
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION IA2IAE(NLXI,NLETA,IA)
      IA2IAE =-1 
      IF (IA.EQ.1) THEN
         IA2IAE = 1 
      ELSEIF (IA.EQ.2) THEN
         IA2IAE = NLXI
      ELSEIF (IA.EQ.3) THEN
         IA2IAE = NLXI*NLETA
      ELSEIF (IA.EQ.4) THEN
         IA2IAE = NLXI*(NLETA - 1) + 1 
      ENDIF
      RETURN
      END 

