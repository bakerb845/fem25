      MODULE INVPLOT_DAT
         INTEGER*4, ALLOCATABLE :: IENGNOD(:) !number of anchor nodes
                                              !on an element 
         INTEGER*4, ALLOCATABLE :: IENGV(:) !vector pointer for IENG 
                                            !in interior domain
         INTEGER*4, ALLOCATABLE :: MASK(:)  !mask(inpg) = inpg_ii or 0
         INTEGER*4 NELEM_II   !number of interior elements
         INTEGER*4 NIENGV     !number of elements in iengv
         INTEGER*4 NNPG_II    !number of anchor nodes in interior
      END MODULE INVPLOT_DAT

      SUBROUTINE DUMP_XZINV_DAT(PROJNM,MGNOD,NNPG, NELEM, NA35,
     ;                          NVINV,NGNOD,LSURF, 
     ;                          CDOMAIN,MASKG,IENG,
     ;                          XLOCS,ZLOCS, WMASK,GRAD)
!
!     Dumps gradient and location pertinent for jacobian regularization
!     search
      USE INVPLOT_DAT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: PROJNM
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG)
      REAL*4, INTENT(IN) :: GRAD(NA35), WMASK(NNPG)
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MASKG(NNPG), 
     ;                         MGNOD, NNPG, NELEM, NVINV, NGNOD, NA35
      LOGICAL*4, INTENT(IN) :: LSURF
!.... local variables
      CHARACTER(80) FLNAME
      REAL*4, ALLOCATABLE :: DAT(:), DAT2(:)
      INTEGER*4 IVEC(4), IA, INPG, INPG_II, IELEM, INPINV, LOC2, IVINV,
     ;          LOC1
!
!----------------------------------------------------------------------!
!
      IF (LSURF) THEN
         FLNAME = TRIM(ADJUSTL(PROJNM))//'_surf_xzinv.dat'
      ELSE
         FLNAME = TRIM(ADJUSTL(PROJNM))//'_body_xzinv.dat'
      ENDIF
      OPEN(UNIT=64,FILE=TRIM(FLNAME),STATUS='REPLACE')
      CALL MASK_IENG(MGNOD,NNPG,NELEM,NGNOD, CDOMAIN,IENG)
      ALLOCATE(DAT(NVINV))
      ALLOCATE(DAT2(NVINV))
      WRITE(64,*) NA35,NNPG_II,NELEM_II,NGNOD,NVINV
      LOC1 = 0
      DO 1 IELEM=1,NELEM_II
         IVEC(1:4) = 0
         DO 2 IA=1,IENGNOD(IELEM)
            LOC1 = LOC1 + 1
            IVEC(IA) = IENGV(LOC1)
    2    CONTINUE
         WRITE(64,840) IENGNOD(IELEM), IVEC(1:NGNOD)
  840    FORMAT(I3,4I10)
    1 CONTINUE
! 
!.... copy over nodal locations and gradient
      DO 3 INPG=1,NNPG
         INPG_II = MASK(INPG)
         IF (INPG_II.GT.0) THEN
            INPINV = MASKG(INPG)
            IF (INPINV.GT.0) THEN
               DO 4 IVINV=1,NVINV
                  LOC2 = (INPINV - 1)*NVINV + IVINV
                  DAT(IVINV) = GRAD(LOC2)
                  DAT2(IVINV) = WMASK(LOC2)
    4          CONTINUE
            ELSE
               DAT(1:NVINV) = 0.0 
               DAT2(1:NVINV) = 0.0
            ENDIF
            IF (NVINV == 1) THEN
               WRITE(64,901) INPG_II,INPINV, 
     ;                     SNGL(XLOCS(INPG)),SNGL(ZLOCS(INPG)),
     ;                     DAT(1:NVINV), DAT2(1:NVINV) 
            ELSEIF (NVINV == 2) THEN
               WRITE(64,902) INPG_II,INPINV, 
     ;                     SNGL(XLOCS(INPG)),SNGL(ZLOCS(INPG)),
     ;                     DAT(1:NVINV), DAT2(1:NVINV) 
            ELSE
               WRITE(64,903) INPG_II,INPINV, 
     ;                     SNGL(XLOCS(INPG)),SNGL(ZLOCS(INPG)),
     ;                     DAT(1:NVINV), DAT2(1:NVINV) 
            ENDIF 
         ENDIF
    3 CONTINUE
  901 FORMAT(I7,I7,2F14.4,E14.6,F12.4)
  902 FORMAT(I7,I7,2F14.4,2E14.6,2F12.4)
  903 FORMAT(I7,I7,2F14.4,3E14.6,3F12.4)
      DEALLOCATE(MASK)
      DEALLOCATE(IENGV)
      DEALLOCATE(IENGNOD)
      DEALLOCATE(DAT) 
      DEALLOCATE(DAT2)
      CLOSE(64)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DUMP_OE_EST(PROJNM, MDIM,MFREQ,MREC, 
     ;                       NDIM,NFREQ,NREC,NSRC, 
     ;                       NFREQ_SRF_INV,NFREQ_BDY_INV, 
     ;                       NSRC_SRF,NSRC_BDY, IRESTP, 
     ;                       AZMOD, PYTAB, CFTYPE,SRCTYP,
     ;                       FREQ, WGHTS,OBS,EST) 
!
!     Dumps data variables pertinent to jacobian regularization search 
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: PROJNM
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*), 
     ;                         EST(MDIM,MFREQ,MREC,*)
      REAL*8, INTENT(IN) :: FREQ(NFREQ), PYTAB(MFREQ,NSRC), AZMOD
      REAL*4, INTENT(IN) :: WGHTS(MDIM,MFREQ,MREC,*)
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                         NFREQ_SRF_INV, NFREQ_BDY_INV, IRESTP, 
     ;                         NSRC_SRF, NSRC_BDY
!.... local variables
      CHARACTER(80) FLNAME
      COMPLEX*8 U, V, W, CN, CE, CZ
      REAL*4 WGT
      INTEGER*4 IFREQ, ISRC, IREC, I, IUNIT
      PARAMETER(IUNIT = 64)
!
!----------------------------------------------------------------------!
!
      FLNAME(1:80) = ' '
      FLNAME = TRIM(ADJUSTL(PROJNM))//'_resid.dat'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      WRITE(IUNIT,*) NFREQ_SRF_INV, NFREQ_BDY_INV
      WRITE(IUNIT,*) NSRC_SRF, NSRC_BDY
      WRITE(IUNIT,*) NREC, IRESTP
      DO 1 IFREQ=1,NFREQ
         WRITE(IUNIT,800) CFTYPE(IFREQ)
  800    FORMAT(A1)
    1 CONTINUE
      DO 2 ISRC=1,NSRC
         WRITE(IUNIT,810) SRCTYP(ISRC)
  810    FORMAT(A2)
    2 CONTINUE
      DO 3 IFREQ=1,NFREQ
         DO 4 ISRC=1,NSRC
            WRITE(IUNIT,*) PYTAB(IFREQ,ISRC)
    4    CONTINUE 
    3 CONTINUE
      DO 5 IFREQ=1,NFREQ
         WRITE(IUNIT,*) FREQ(IFREQ)
         DO 6 ISRC=1,NSRC
            DO 7 IREC=1,NREC
               U = EST(1,IFREQ,IREC,ISRC)
               V = EST(2,IFREQ,IREC,ISRC)
               W = EST(3,IFREQ,IREC,ISRC)
               CALL CROTATE(+SNGL(AZMOD),U,V, CN,CE)
               CZ = W
               DO 8 I=1,NDIM
                  WGT = WGHTS(I,IFREQ,IREC,ISRC)
                  WRITE(IUNIT,900) 
     ;            OBS(I,IFREQ,IREC,ISRC),wgt,i,ifreq,irec,isrc,'o'
                  IF (I.EQ.1) THEN
                     IF (WGT > 0.0) THEN
                        WRITE(IUNIT,900) U,1.0,i,ifreq,irec,isrc,'e' 
                     ELSE
                        WRITE(IUNIT,900) U,0.0,i,ifreq,irec,isrc,'e'
                     ENDIF
                  ELSEIF (I.EQ.2) THEN
                     IF (WGT > 0.0) THEN
                        WRITE(IUNIT,900) V,1.0,i,ifreq,irec,isrc,'e'
                     ELSE
                        WRITE(IUNIT,900) V,0.0,i,ifreq,irec,isrc,'e'
                     ENDIF
                  ELSE
                     IF (WGT > 0.0) THEN
                        WRITE(IUNIT,900) W,1.0,i,ifreq,irec,isrc,'e'
                     ELSE
                        WRITE(IUNIT,900) W,0.0,i,ifreq,irec,isrc,'e'
                     ENDIF
                  ENDIF
    8          CONTINUE
    7       CONTINUE
    6    CONTINUE
    5 CONTINUE
      CLOSE(IUNIT)
  900 FORMAT('(',E14.5,',',E14.5,')',F12.5,I2,I4,I4,I4,2X,A1)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_SHGRAD_VTK(PROJNM,MGNOD,NNPG, NELEM, 
     ;                           CINVTYPE, ITYPE, IBLOCK,K,IALPHA, 
     ;                           CDOMAIN,MASKG,IENG, 
     ;                           XLOCS,ZLOCS, SHGRAD)
! 
!     Plots the gradients during the inversion.  The output is a 
!     binary .vtk file.  Naming scheme is 
!       ./invfiles_vtk/projnm-block-iteration-step_g??.vtk 
! 
!     INPUT      MEANING
!     -----      ------- 
!     CDOMAIN    'I' interior; 'E' bielak; 'A' absorbing layer
!     CINVTYPE   inversion type 'PP', 'SS', or 'PS' right now
!     IALPHA     step length counter
!     IBLOCK     block number
!     IENG       connectivity vector for anchor nodes
!     ITYPE      = 1 Gradient, = 2 Hessian, = 3 search direction, 
!                = 4 approximate newton direction, 5 gauss-newton
!     K          iteration in iterative loop
!     MGNOD      max number of anchor nodes
!     NELEM      number of elements 
!     NNPG       number of anchor nodes 
!     PROJNM     project name
!     SHGRAD     search, Hessian, or gradient direction
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!    
!.... variable declarations
      USE INVPLOT_DAT
      CHARACTER(*), INTENT(IN) :: PROJNM
      CHARACTER(2), INTENT(IN) :: CINVTYPE
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG)
      REAL*4, INTENT(IN) :: SHGRAD(*)  
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MASKG(NNPG), MGNOD,
     ;                         NNPG, NELEM,ITYPE, 
     ;                         IBLOCK,K,IALPHA
!.... local variables
      CHARACTER(80) FILENM
      CHARACTER(9) APP 
      CHARACTER(3) CBLOCK,CK,CALPHA 
      CHARACTER(1), ALLOCATABLE :: FPASS(:)
      REAL*4, ALLOCATABLE :: XPTS(:), ZPTS(:), DAT(:) 
c     INTEGER*4, ALLOCATABLE :: IENGV(:), IENGNOD(:)
      LOGICAL*4 LEX, LISDIR 
! 
!----------------------------------------------------------------------!
! 
!.... initial check on cinvtype
      ITYPEP = 0
      APP(1:9) = ' '
      IF (ITYPE.EQ.1) THEN !gradient
         IF (CINVTYPE.EQ.'PP' .OR. CINVTYPE.EQ.'pp') THEN
            APP = '_gpp.vtk'
            NVINV = 1
            ITYPEP = 1
         ELSEIF (CINVTYPE.EQ.'SS' .OR. CINVTYPE.EQ.'ss') THEN
            APP = '_gss.vtk'
            NVINV = 1
            ITYPEP = 2
         ELSEIF (CINVTYPE.EQ.'PS' .OR. CINVTYPE.EQ.'ps') THEN
            APP = '_gps.vtk'
            NVINV = 2
            ITYPEP = 1
         ELSE
            WRITE(*,*) 'plot_shgrad_vtk: Cant classify gradient type!'
            RETURN
         ENDIF
      ELSEIF(ITYPE.EQ.2) THEN !hessian
         IF (CINVTYPE.EQ.'PP' .OR. CINVTYPE.EQ.'pp') THEN
            APP = '_hpp.vtk'
            NVINV = 1
            ITYPEP = 3
         ELSEIF (CINVTYPE.EQ.'SS' .OR. CINVTYPE.EQ.'ss') THEN
            APP = '_hss.vtk'
            NVINV = 1
            ITYPEP = 4
         ELSEIF (CINVTYPE.EQ.'PS' .OR. CINVTYPE.EQ.'ps') THEN
            APP = '_hps.vtk'
            NVINV = 2
            ITYPEP = 2
         ELSEIF (CINVTYPE.EQ.'SP' .OR. CINVTYPE.EQ.'sp') THEN
            APP = '_hsp.vtk'
            NVINV = 1
            ITYPEP = 7 
         ELSE
            WRITE(*,*) 'plot_shgrad_vtk: Cant classify gradient type!'
            RETURN
         ENDIF
      ELSEIF (ITYPE.EQ.3) THEN !search direction 
         IF (CINVTYPE.EQ.'PP' .OR. CINVTYPE.EQ.'pp') THEN
            APP = '_spp.vtk'
            NVINV = 1 
            ITYPEP = 5
         ELSEIF (CINVTYPE.EQ.'SS' .OR. CINVTYPE.EQ.'ss') THEN
            APP = '_sss.vtk'
            NVINV = 1 
            ITYPEP = 6
         ELSEIF (CINVTYPE.EQ.'PS' .OR. CINVTYPE.EQ.'ps') THEN
            APP = '_sps.vtk'
            NVINV = 2 
            ITYPEP = 3
         ELSE
            WRITE(*,*) 'plot_shgrad_vtk: Cant classify gradient type!'
            RETURN
         ENDIF
      ELSEIF (ITYPE.EQ.4) THEN !higher order scatteringnewton direction
         IF (CINVTYPE.EQ.'PP' .OR. CINVTYPE.EQ.'pp') THEN
            APP = '_p1pp.vtk'
            NVINV = 1
            ITYPEP = 5
         ELSEIF (CINVTYPE.EQ.'SS' .OR. CINVTYPE.EQ.'ss') THEN
            APP = '_p1ss.vtk'
            NVINV = 1
            ITYPEP = 6
         ELSEIF (CINVTYPE.EQ.'PS' .OR. CINVTYPE.EQ.'ps') THEN
            APP = '_p1ps.vtk'
            NVINV = 2
            ITYPEP = 3
         ELSE
            WRITE(*,*) 'plot_shgrad_vtk: Cant classify gradient type!'
            RETURN
         ENDIF
      ELSE !IF (ITYPE.EQ.5) THEN !higher order scatteringnewton direction
         IF (CINVTYPE.EQ.'PP' .OR. CINVTYPE.EQ.'pp') THEN
            APP = '_p2pp.vtk'
            NVINV = 1 
            ITYPEP = 5 
         ELSEIF (CINVTYPE.EQ.'SS' .OR. CINVTYPE.EQ.'ss') THEN
            APP = '_p2ss.vtk'
            NVINV = 1 
            ITYPEP = 6 
         ELSEIF (CINVTYPE.EQ.'PS' .OR. CINVTYPE.EQ.'ps') THEN
            APP = '_p2ps.vtk'
            NVINV = 2 
            ITYPEP = 3 
         ELSE
            WRITE(*,*) 'plot_shgrad_vtk: Cant classify gradient type!'
            RETURN
         ENDIF
      ENDIF
      APP = ADJUSTL(APP) 
!.... extract mesh partition and convert IENG to a vector 
c     NWORK = 4*NELEM          !worst case is all linear quads 
c     ALLOCATE(IENGV(NWORK))   !will hold vectorized IENG vector 
c     ALLOCATE(IENGNOD(NELEM)) !eventually will be input
      CALL MASK_IENG(MGNOD,NNPG,NELEM,4, CDOMAIN,IENG)
      ALLOCATE(XPTS(NNPG_II))     !VisIt uses single precision
      ALLOCATE(ZPTS(NNPG_II))     !VisIt uses single precision
! 
!.... copy over nodal locations and gradient
      ALLOCATE(DAT(NVINV*NNPG_II))
      DAT(1:NVINV*NNPG_II) = 0.0
      DO 1 INPG=1,NNPG
         INPG_II = MASK(INPG)
         IF (INPG_II.GT.0) THEN
            XPTS(INPG_II) = SNGL(XLOCS(INPG))
            ZPTS(INPG_II) = SNGL(ZLOCS(INPG))
            INPINV = MASKG(INPG)
            IF (INPINV.GT.0) THEN
               DO 2 IVINV=1,NVINV
                  LOC1 = (IVINV - 1)*NNPG_II + INPG_II
                  !LOC2 = (IVINV - 1)*NNPINV + INPINV
                  LOC2 = (INPINV - 1)*NVINV + IVINV
                  DAT(LOC1) = SHGRAD(LOC2)
    2          CONTINUE
            ENDIF
         ENDIF
    1 CONTINUE
!
!.... check if directory exists
      LEX = LISDIR('./invfiles_vtk')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./invfiles_vtk')
!
!.... file name
      CBLOCK(1:3) = ' '
      CK    (1:3) = ' '
      CALPHA(1:3) = ' '
      WRITE(CBLOCK,'(I3)') IBLOCK
      WRITE(CK    ,'(I3)') K
      WRITE(CALPHA,'(I3)') IALPHA
      CBLOCK = ADJUSTL(CBLOCK)
      CK = ADJUSTL(CK)
      CALPHA = ADJUSTL(CALPHA)
      CALL NULLS(80,FILENM)
!
!.... step length run
      IF (IALPHA.GT.0) THEN
         FILENM = './invfiles_vtk/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//
     ;            '-'//TRIM(CK)//'-'//CALPHA//TRIM(APP)
      ELSE
         FILENM = './invfiles_vtk/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//
     ;            '-'//TRIM(CK)//TRIM(APP)
      ENDIF
      FILENM = ADJUSTL(FILENM)
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 5 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    5 CONTINUE
      IF (NVINV.EQ.1) THEN
         CALL VTK_SHGRAD1(FPASS,LENFL, NNPG_II,NNPG_II,NELEM_II,NIENGV,
     ;                    ITYPEP,IENGNOD,IENGV, XPTS,ZPTS,DAT)
      ELSEIF (NVINV.EQ.2) THEN
         CALL VTK_SHGRAD2(FPASS,LENFL, NNPG_II,NNPG_II,NELEM_II,NIENGV,
     ;                    ITYPEP,IENGNOD,IENGV, XPTS,ZPTS,DAT)
      ELSE
         WRITE(*,*) 'plot_hgrad_vtk: I dont know how to write nvinv > 2'
      ENDIF
      DEALLOCATE(MASK)
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
      SUBROUTINE PLOT_ELMOD_UPD(PROJNM,MGNOD,MNPG,NNPG, NELEM, LISISO, 
     ;                          IBLOCK,K,IALPHA, CDOMAIN, 
     ;                          IENG,XLOCS,ZLOCS, DENS,ECOEFF)
! 
!     Plots the elastic model.  The output is a binary .vtk file 
! 
!     INPUT      MEANING
!     -----      ------- 
!     DENS       density at anchor nodes
!     ECOEFF     elastic coefficients at anchor nodes
!     IALPHA     step length iteration number (0 beginning steplen)
!     IBLOCK     inversion block number
!     IENG       connectivity vector for anchor nodes
!     K          iteraiton number
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
      USE INVPLOT_DAT
      CHARACTER(*), INTENT(IN) :: PROJNM
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: ECOEFF(MNPG,*), DENS(NNPG), 
     ;                      XLOCS(NNPG), ZLOCS(NNPG)
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MGNOD,MNPG,NNPG, 
     ;                         NELEM, IBLOCK,K,IALPHA 
      LOGICAL*4, INTENT(IN) :: LISISO 
!.... local variables
      CHARACTER(4), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FILENM
      CHARACTER(3) CBLOCK,CK,CALPHA
      CHARACTER(1), ALLOCATABLE :: FPASS(:)
      REAL*8 ALPHA2, BETA2 
      REAL*4, ALLOCATABLE :: XPTS(:), ZPTS(:), DAT(:) 
c     INTEGER*4, ALLOCATABLE :: IENGV(:), IENGNOD(:)
      LOGICAL*4 LEX, LISDIR, LSWAP  
      CHARACTER(4) PACKR4, PACKI4  
      INTEGER*4 NGNOD4, ENDIAN
      PARAMETER(IUNIT = 65)
! 
!----------------------------------------------------------------------!
! 
!.... extract mesh partition and convert IENG to a vector 
      NWORK = 4*NELEM          !worst case is all linear quads 
      !ALLOCATE(IENGV(NWORK))   !will hold vectorized IENG vector 
      !ALLOCATE(IENGNOD(NELEM)) !eventually will be input
      CALL MASK_IENG(MGNOD,NNPG,NELEM,4, CDOMAIN,IENG)
      ALLOCATE(XPTS(NNPG_II))     !VisIt uses single precision
      ALLOCATE(ZPTS(NNPG_II))     !VisIt uses single precision
! 
!.... copy over nodal locations
      IF (LISISO) THEN
         ALLOCATE(DAT(3*NNPG))
      ELSE
         WRITE(*,*) 'plot_model_vtk: Error program lisiso = false'
      ENDIF
      DO 1 INPG=1,NNPG
         INPG_II = MASK(INPG)
         IF (INPG_II.GT.0) THEN
            XPTS(INPG_II) = SNGL(XLOCS(INPG))
            ZPTS(INPG_II) = SNGL(ZLOCS(INPG))
            IF (LISISO) THEN
               LOC1 =             INPG_II
               LOC2 =   NNPG_II + INPG_II
               LOC3 = 2*NNPG_II + INPG_II
               ALPHA2 = (ECOEFF(INPG,1)+2.D0*ECOEFF(INPG,2))/DENS(INPG)
               BETA2 = ECOEFF(INPG,2)/DENS(INPG)
               DAT(LOC1) = SNGL(DSQRT(ALPHA2))
               DAT(LOC2) = SNGL(DSQRT(BETA2))
               DAT(LOC3) = SNGL(DENS(INPG))
            ELSE
               WRITE(*,*) 'plot_model_vtk: Error no anistropic ordering'
               RETURN
            ENDIF
         ENDIF !end check on active node
    1 CONTINUE
! 
!.... write vtk file
      LEX = LISDIR('./invmodel_vtk')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./invmodel_vtk')
!
!.... file name
      CBLOCK(1:3) = ' ' 
      CK    (1:3) = ' ' 
      CALPHA(1:3) = ' ' 
      WRITE(CBLOCK,'(I3)') IBLOCK
      WRITE(CK    ,'(I3)') K
      WRITE(CALPHA,'(I3)') IALPHA
      CBLOCK = ADJUSTL(CBLOCK)
      CK = ADJUSTL(CK)
      CALPHA = ADJUSTL(CALPHA)
      CALL NULLS(80,FILENM)
      IF (IALPHA.GT.0) THEN
         FILENM = './invmodel_vtk/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//
     ;            '-'//TRIM(CK)//'-'//CALPHA//'_invmod.vtk'
      ELSE
         FILENM = './invmodel_vtk/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//
     ;            '-'//TRIM(CK)//'_invmod.vtk'
      ENDIF
      FILENM = ADJUSTL(FILENM)
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 2 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    2 CONTINUE
      CALL VTK_MODEL(FPASS,LENFL, NNPG_II,NNPG_II,NELEM_II,NIENGV,
     ;               IENGNOD,IENGV, XPTS,ZPTS,DAT)
!
!.... and dump model
      FILENM(1:80) = ' '
      IF (IALPHA.GT.0) THEN
         FILENM = './invmodel_vtk/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//
     ;            '-'//TRIM(CK)//'-'//TRIM(CALPHA)//'_invmod.bin'
      ELSE
         FILENM = './invmodel_vtk/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//
     ;            '-'//TRIM(CK)//'_invmod.bin'
      ENDIF
      FILENM = ADJUSTL(FILENM)
!
!.... determine endianness 
      MYEND = ENDIAN()
      LSWAP = .FALSE.
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... pack the data
      NGNOD4 = 4
      LWORK = 4 + 2*NELEM_II*NGNOD4 + 5*NNPG_II
      ALLOCATE(CDAT(LWORK))
      CDAT(1) = PACKI4(LSWAP,NNPG_II) 
      CDAT(2) = PACKI4(LSWAP,NELEM_II)
      CDAT(3) = PACKI4(LSWAP,NGNOD4) 
      IF (.NOT.LISISO) THEN
         WRITE(*,* )'Error dont know how to write model file'
         RETURN
      ENDIF
      CDAT(4) = PACKI4(LSWAP,3)
      INDX = 4
      JNDX = 0
      DO 5 IELEM=1,NELEM
         IF (CDOMAIN(IELEM) == 'I') THEN 
            DO 6 IA=1,NGNOD4
               !global anchor node location
               INDX = INDX + 1
               CDAT(INDX) = PACKI4(LSWAP,IENG(IA,IELEM))
               !connect to local pointers
               INDX = INDX + 1
               IF (IENG(IA,IELEM) > 0) THEN !
                  JNDX = JNDX + 1
                  CDAT(INDX) = PACKI4(LSWAP,IENGV(JNDX))
               ELSE
                  CDAT(INDX) = PACKI4(LSWAP,0)
               ENDIF
    6       CONTINUE
         ENDIF
    5 CONTINUE
      DO 7 INPG_II=1,NNPG_II
         INDX = INDX + 1
         CDAT(INDX) = PACKR4(LSWAP,XPTS(INPG_II))
         INDX = INDX + 1
         CDAT(INDX) = PACKR4(LSWAP,ZPTS(INPG_II))
         LOC1 = INPG_II
         LOC2 = NNPG_II + INPG_II
         LOC3 = 2*NNPG_II + INPG_II
         INDX = INDX + 1
         CDAT(INDX) = PACKR4(LSWAP,DAT(LOC1))
         INDX = INDX + 1
         CDAT(INDX) = PACKR4(LSWAP,DAT(LOC2))
         INDX = INDX + 1
         CDAT(INDX) = PACKR4(LSWAP,DAT(LOC3)) 
    7 CONTINUE 
      NBYTES = LWORK*4 
      OPEN(UNIT=IUNIT,FILE=TRIM(FILENM),FORM='UNFORMATTED',
     ;     STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      WRITE(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,LWORK) 
      CLOSE(IUNIT)
      DEALLOCATE(MASK)
      DEALLOCATE(CDAT) 
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
      SUBROUTINE PLOT_JACOB_SB(MDIM,MGNOD, NCOMP,NNPG,NELEM,IREC,ISRC, 
     ;                         LSURF,CDOMAIN,MASKG,IENG, 
     ;                         RJAC,QJAC, FREQ, XLOCS,ZLOCS) 
!
!     Plots the complex valued Jacobian for a frequency, source, 
!     receiver pair
      USE INVPLOT_DAT
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), FREQ
      REAL*4, INTENT(IN) :: RJAC(MDIM,*), QJAC(MDIM,*) 
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MASKG(NNPG),  MDIM, MGNOD,
     ;                         NCOMP, NNPG, NELEM, IREC, ISRC
      LOGICAL*4, INTENT(IN) :: LSURF
!.... local variables
      CHARACTER(1), ALLOCATABLE :: FPASS(:)
      REAL*4, ALLOCATABLE :: XPTS(:), ZPTS(:), RMAG(:), PHASE(:) 
      CHARACTER(80) FILENM
      CHARACTER(12) CFREQ
      CHARACTER(5) CREC, CSRC 
      REAL*4 PI180I 
      INTEGER*4 NWORK, INPG, INPINV, I, LOC1, LOC2, 
     ;          LENFL, LDD, INPG_II 
      LOGICAL*4 LISDIR, LEX 
      PARAMETER(PI180I = 57.29577951308232) 
!
!----------------------------------------------------------------------!

!.... extract mesh partition and convert IENG to a vector 
      NWORK = 4*NELEM          !worst case is all linear quads 
      CALL MASK_IENG(MGNOD,NNPG,NELEM,4, CDOMAIN,IENG)
      LDD   = NNPG_II
      ALLOCATE(XPTS(NNPG_II))   !VisIt uses single precision
      ALLOCATE(ZPTS(NNPG_II))   !VisIt uses single precision
! 
!.... copy over nodal locations and gradient
      ALLOCATE(RMAG(NCOMP*LDD))
      ALLOCATE(PHASE(NCOMP*LDD))
      RMAG(1:NCOMP*LDD) = 0.0
      PHASE(1:NCOMP*LDD) = 0.0
      DO 1 INPG=1,NNPG
         INPG_II = MASK(INPG)
         IF (INPG_II.GT.0) THEN
            XPTS(INPG_II) = SNGL(XLOCS(INPG))
            ZPTS(INPG_II) = SNGL(ZLOCS(INPG))
            INPINV = MASKG(INPG)
            IF (INPINV.GT.0) THEN
               DO 2 I=1,NCOMP
                  LOC1 = (I - 1)*LDD + INPG_II
                  RMAG(LOC1)  = CABS(CMPLX(RJAC(I,INPINV),
     ;                                     QJAC(I,INPINV)))
                  LOC2 = (I - 1)*LDD + INPG_II
                  PHASE(LOC2) = ATAN2(QJAC(I,INPINV),
     ;                                RJAC(I,INPINV))*PI180I
    2          CONTINUE
            ENDIF
         ENDIF
    1 CONTINUE
!
!.... name file and write
      LEX = LISDIR('./jacobian_vtk')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./jacobian_vtk')
      FILENM(1:80) = ' '
      CFREQ(1:12) = ' '
      WRITE(CFREQ,'(F12.5)') FREQ
      CFREQ = ADJUSTL(CFREQ)
      CSRC(1:5) = ' '
      WRITE(CSRC,'(I5)') ISRC
      CSRC = ADJUSTL(CSRC)
      CREC(1:5) = ' '
      WRITE(CREC,'(I5)') IREC
      CREC = ADJUSTL(CREC)
      IF (LSURF) THEN
         FILENM = './jacobian_vtk/jacobian_srf-'//TRIM(CFREQ)//'-'//
     ;            TRIM(CSRC)//'-'//TRIM(CREC)//'.vtk'
      ELSE
         FILENM = './jacobian_vtk/jacobian_bdy-'//TRIM(CFREQ)//'-'//
     ;            TRIM(CSRC)//'-'//TRIM(CREC)//'.vtk'
      ENDIF
      FILENM = ADJUSTL(FILENM)
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 3 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    3 CONTINUE
!
!.... write file
      CALL VTK_QRESP25D_EL(FPASS,LENFL, LDD,NNPG_II,NELEM_II,NIENGV,
     ;                     IENGNOD,IENGV, XPTS,ZPTS,RMAG,PHASE)
!
!.... clean space
      DEALLOCATE(MASK)
      DEALLOCATE(IENGV)
      DEALLOCATE(IENGNOD)
      DEALLOCATE(XPTS)
      DEALLOCATE(ZPTS)
      DEALLOCATE(RMAG)
      DEALLOCATE(PHASE)
      DEALLOCATE(FPASS)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_JACOB(MDIM,MGNOD, NCOMP,NNPG,NELEM,IREC,ISRC, 
     ;                      CDOMAIN,MASKG,IENG, RJAC,QJAC,  
     ;                      FREQ, XLOCS,ZLOCS) 
!
!     Plots the complex valued Jacobian for a frequency, source, 
!     receiver pair

      USE INVPLOT_DAT
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), FREQ
      REAL*4, INTENT(IN) :: RJAC(MDIM,*), QJAC(MDIM,*) 
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MASKG(NNPG),  MDIM, MGNOD,
     ;                         NCOMP, NNPG, NELEM, IREC, ISRC
!.... local variables
      CHARACTER(1), ALLOCATABLE :: FPASS(:)
      REAL*4, ALLOCATABLE :: XPTS(:), ZPTS(:), RMAG(:), PHASE(:) 
c     INTEGER*4, ALLOCATABLE :: IENGV(:), IENGNOD(:) 
      CHARACTER(80) FILENM
      CHARACTER(12) CFREQ
      CHARACTER(5) CREC, CSRC 
      REAL*4 PI180I 
      INTEGER*4 NWORK, INPG, INPINV, I, LOC1, LENFL, LDD, INPG_II 
      LOGICAL*4 LISDIR, LEX 
      PARAMETER(PI180I = 57.29577951308232) 
!
!----------------------------------------------------------------------!

!.... extract mesh partition and convert IENG to a vector 
      NWORK = 4*NELEM          !worst case is all linear quads 
      CALL MASK_IENG(MGNOD,NNPG,NELEM,4, CDOMAIN,IENG)
      LDD   = NNPG_II
      ALLOCATE(XPTS(NNPG_II))   !VisIt uses single precision
      ALLOCATE(ZPTS(NNPG_II))   !VisIt uses single precision
! 
!.... copy over nodal locations and gradient
      ALLOCATE(RMAG(NCOMP*LDD)) 
      ALLOCATE(PHASE(NCOMP*LDD))
      RMAG(1:NCOMP*LDD) = 0.0
      PHASE(1:NCOMP*LDD) = 0.0 
      DO 1 INPG=1,NNPG
         INPG_II = MASK(INPG)
         IF (INPG_II.GT.0) THEN
            XPTS(INPG_II) = SNGL(XLOCS(INPG))
            ZPTS(INPG_II) = SNGL(ZLOCS(INPG))
            INPINV = MASKG(INPG)
            IF (INPINV.GT.0) THEN
               DO 2 I=1,NCOMP 
                  LOC1 = (I - 1)*LDD + INPG_II
                  RMAG(LOC1)  = CABS(CMPLX(RJAC(I,INPINV),
     ;                                     QJAC(I,INPINV)))
                  PHASE(LOC1) = ATAN2(QJAC(I,INPINV),
     ;                                RJAC(I,INPINV))*PI180I 
    2          CONTINUE    
            ENDIF
         ENDIF
    1 CONTINUE

      LEX = LISDIR('./jacobian_vtk')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./jacobian_vtk')
      FILENM(1:80) = ' ' 
      CFREQ(1:12) = ' '
      WRITE(CFREQ,'(F12.5)') FREQ
      CFREQ = ADJUSTL(CFREQ) 
      CSRC(1:5) = ' '
      WRITE(CSRC,'(I5)') ISRC
      CSRC = ADJUSTL(CSRC)
      CREC(1:5) = ' '
      WRITE(CREC,'(I5)') IREC 
      CREC = ADJUSTL(CREC)
      FILENM = './jacobian_vtk/jacobian-'//TRIM(CFREQ)//'-'//
     ;         TRIM(CSRC)//'-'//TRIM(CREC)//'.vtk'
      FILENM = ADJUSTL(FILENM)
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 3 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    3 CONTINUE
!
!.... write file
      CALL VTK_QRESP25D_EL(FPASS,LENFL, LDD,NNPG_II,NELEM_II,NIENGV,   
     ;                     IENGNOD,IENGV, XPTS,ZPTS,RMAG,PHASE)   
!
!.... clean space
      DEALLOCATE(MASK)
      DEALLOCATE(IENGV)
      DEALLOCATE(IENGNOD)
      DEALLOCATE(XPTS)
      DEALLOCATE(ZPTS)
      DEALLOCATE(RMAG)
      DEALLOCATE(PHASE) 
      DEALLOCATE(FPASS)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_ADJGRAD_VTK(MGNOD,MDIM,MEN, NDIM,NNPG,NDOF,
     ;                            NLXI,NLETA, NELEM,NCOMP, ISRC,ITYPE, 
     ;                            FREQ, LM,IENG, XLOCS,ZLOCS, U)

      COMPLEX*8, INTENT(IN) :: U(NDOF) 
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), FREQ  
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), ITYPE
!.... local variables
      CHARACTER(1), ALLOCATABLE :: FPASS(:) 
      REAL*4, ALLOCATABLE :: RMAG(:), PHASE(:), XPTS(:), ZPTS(:) 
      INTEGER*4, ALLOCATABLE :: IDPTR(:,:), IENGNOD(:), IENGV(:) 
      CHARACTER(80) FILENM 
      CHARACTER(12) CFREQ
      CHARACTER(8) CAPP
      CHARACTER(5) CSRC
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
               PHASE(INDX) = ATAN2(IMAG(U(IDOF)),REAL(U(IDOF)))/PI180
            ENDIF
    4    CONTINUE
    3 CONTINUE
      DEALLOCATE(IDPTR)
      CALL NULLS(80,FILENM)
      CALL NULLS(12,CFREQ)
      CALL NULLS(5,CSRC)
      CALL NULLS(8,CAPP) 
      IF (ITYPE.EQ.1) THEN
         CAPP = '_lam.vtk'
      ELSEIF (ITYPE.EQ.2) THEN
         CAPP = '_nu.vtk'
      ELSEIF (ITYPE.EQ.3) THEN
         CAPP = '_alp.vtk'
      ELSE
         CAPP = '.vtk'
      ENDIF
      CAPP = ADJUSTL(CAPP)  
      WRITE(CFREQ,'(F12.5)') FREQ
      CFREQ = ADJUSTL(CFREQ)
      WRITE(CSRC,'(I5)') ISRC
      CSRC = ADJUSTL(CSRC)
      LEX = LISDIR('./back_prop')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./back_prop')
      FILENM = './back_prop/bpwave-'//
     ;         TRIM(CFREQ)//'-'//TRIM(CSRC)//TRIM(CAPP) !'.vtk'
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
      SUBROUTINE GETBDIAG_HESS(NHSIZE,NNPINV,NVINV,NA35, JOB, 
     ;                         IPIV,GRADPC, HESS) 
!
!     Gets the block diaognal from the Hessian.  Note the Hessian 
!     has been LU factored so we need to extract it 
      REAL*4, INTENT(IN) :: GRADPC(NHSIZE)  
      INTEGER*4, INTENT(IN) :: IPIV(NA35), NHSIZE, NNPINV, NVINV, NA35,
     ;                         JOB
      REAL*4, INTENT(OUT) :: HESS(*) 
      REAL*4, ALLOCATABLE :: BLOCK(:,:), LMAT(:,:), UMAT(:,:) 
      INTEGER*4 IVINV,JVINV,KVINV, INDX,IDEST,INPINV  
!
!----------------------------------------------------------------------!
!
!.... intialize
      ALLOCATE(LMAT(NVINV,NVINV))
      ALLOCATE(UMAT(NVINV,NVINV)) 
      ALLOCATE(BLOCK(NVINV,NVINV)) 
      HESS(1:NNPINV) = 0.0  
      UMAT(1:NVINV,1:NVINV) = 0.0
      LMAT(1:NVINV,1:NVINV) = 0.0
      DO 1 IVINV=1,NVINV
         LMAT(IVINV,IVINV) = 1.D0
    1 CONTINUE 
!
!.... pull out LU blocks -> multiply, extract diagonal
      INDX = 0 
      IDEST = 0
      DO 2 INPINV=1,NNPINV 
         DO 3 KVINV=1,NVINV
            IVINV = IPIV(KVINV)
            DO 4 JVINV=1,NVINV
               INDX = INDX + 1
               IF (IVINV.LT.JVINV) THEN
                  LMAT(IVINV,JVINV) = GRADPC(INDX) 
               ELSE
                  UMAT(IVINV,JVINV) = GRADPC(INDX) 
               ENDIF
    4       CONTINUE
    3    CONTINUE 
         BLOCK = MATMUL(LMAT,UMAT) 
         DO 5 IVINV=1,NVINV
            DO 6 JVINV=1,NVINV  
               IF (JOB.EQ.1 .AND. IVINV.EQ.1 .AND. JVINV.EQ.1) THEN
                  IDEST = IDEST + 1
                  HESS(IDEST) = BLOCK(IVINV,JVINV)
               ELSEIF (JOB.EQ.1 .AND. IVINV.EQ.2 .AND. JVINV.EQ.2) THEN
                  IDEST = IDEST + 1
                  HESS(IDEST) = BLOCK(IVINV,JVINV)
               ELSEIF (JOB.EQ.2 .AND. IVINV.EQ.1 .AND. JVINV.EQ.2) THEN
                  IDEST = IDEST + 1
                  HESS(IDEST) = BLOCK(IVINV,JVINV) 
               ENDIF  
    6       CONTINUE
    5    CONTINUE
    2 CONTINUE !loop on nodal points in inversion
!
!.... clean space
      IF (ALLOCATED(UMAT)) DEALLOCATE(UMAT)
      IF (ALLOCATED(LMAT)) DEALLOCATE(LMAT)
      IF (ALLOCATED(BLOCK)) DEALLOCATE(BLOCK)           
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PLOT_APHESS_VTK(PROJNM, LDH,N, IBLOCK,K, HESS)
!
!     Since Gnuplot wants to self destruct while plotting the 
!     approximate Hessian, Re( adj(J) J), I will let vtk do it
! 
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM 
      REAL*4, INTENT(IN) :: HESS(LDH,*) 
      INTEGER*4, INTENT(IN) :: LDH,N, IBLOCK, K  
!.... local variables
      CHARACTER(80) FILENM
      CHARACTER(3) CBLOCK, CK
      CHARACTER(1), ALLOCATABLE :: FPASS(:)
      REAL*4, ALLOCATABLE :: HESSB(:) 
      LOGICAL*4 LISDIR, LEX 
!
!----------------------------------------------------------------------!
!
!.... set space and copy over buffer for C
      NWORK = N**2
      ALLOCATE(HESSB(NWORK),STAT=IERR) 
      IF (IERR.NE.0) THEN
         WRITE(*,*) 'plot_aphess_vtk: Error setting space!'
         RETURN
      ENDIF
      INDX = 0 
      DO 1 I=1,N
         DO 2 J=N,1,-1 !flips up and down 
            INDX = INDX + 1
            HESSB(INDX) = HESS(I,J) 
    2    CONTINUE 
    1 CONTINUE 
!
!.... set file name
      LEX = LISDIR('./invfiles_vtk')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./invfiles_vtk') 
      FILENM(1:80) = ' '
      CBLOCK(1:3) = ' '
      CK(1:3) = ' '
      WRITE(CBLOCK,'(I3)') IBLOCK 
      WRITE(CK,'(I3)') K 
      FILENM = './invfiles_vtk/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//
     ;         '-'//TRIM(CK)//'_hessian.vtk'
      FILENM = ADJUSTL(FILENM)
      CALL DBLANK(80,FILENM,LENFL)
      ALLOCATE(FPASS(LENFL))
      DO 4 I=1,LENFL
         FPASS(I) = FILENM(I:I)
    4 CONTINUE
!
!.... and some parameters for VisIt
      NX = N
      NY = 1 !Hessian is a matrix, so this is a 1D thing
      NZ = N
      CALL VTK_PLOT_REG_HESS(FPASS,LENFL, N, HESSB) 
      DEALLOCATE(HESSB) 
      DEALLOCATE(FPASS) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE MASK_IENG(MGNOD,NNPG,NELEM,NGNOD, CDOMAIN,IENG)
!
!     When plotting the inverse problem we don't really care what 
!     happens in the bielak layer or absorbing layer, so we truncated
!     the domain 
!
!.... variable declarations
      USE INVPLOT_DAT
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MGNOD,NNPG,NELEM, NGNOD
!.... local variables
      INTEGER*4 IELEM, IA, INPG, INPG_II
      LOGICAL*4, ALLOCATABLE :: LACTIVE(:) 
!
!----------------------------------------------------------------------!
!
      ALLOCATE(LACTIVE(NNPG))
      NIENGV = 0
      NELEM_II = 0
      LACTIVE(1:NNPG) = .FALSE.
      DO 1 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'I') THEN
            NELEM_II = NELEM_II + 1
            DO 2 IA=1,NGNOD 
               INPG = IENG(IA,IELEM)
               IF (INPG.EQ.0) GOTO 20
               NIENGV = NIENGV + 1
               LACTIVE(INPG) = .TRUE.
    2       CONTINUE
   20       CONTINUE 
         ENDIF
    1 CONTINUE 
!
!.... for the active elements save IENGV
      NNPG_II = 0
      ALLOCATE(MASK(NNPG))
      MASK(1:NNPG) = 0
      DO 3 INPG=1,NNPG 
         IF (LACTIVE(INPG)) THEN
            NNPG_II = NNPG_II + 1
            MASK(INPG) = NNPG_II
         ENDIF
    3 CONTINUE 
!
!.... generate the IENGV vector
      ALLOCATE(IENGV(NIENGV))
      ALLOCATE(IENGNOD(NELEM_II))
      IELEM_II = 0
      NIENGV = 0
      DO 4 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'I') THEN
            IELEM_II = IELEM_II + 1
            IENGNOD(IELEM_II) = 4
            DO 5 IA=1,NGNOD
               INPG = IENG(IA,IELEM)
               IF (INPG.EQ.0) THEN
                  IENGNOD(IELEM_II) = 3
                  GOTO 50
               ENDIF
               NIENGV = NIENGV + 1
               INPG_II = MASK(INPG)
               IENGV(NIENGV) = INPG_II
    5       CONTINUE
   50       CONTINUE 
         ENDIF 
    4 CONTINUE  
      DEALLOCATE(LACTIVE)
      RETURN
      END 
