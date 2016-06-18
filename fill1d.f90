      SUBROUTINE FILL1D(MSH,MOD1D)
!
!     Fills in the headers on the 1D models.  Just saves space in the 
!     driver programs. 
!
!     INPUT      MEANING
!     -----      ------- 
!     MSH        data structure with mesh parameters and 2D model values
!
!     OUTPUT     MEANING
!     ------     ------- 
!     HDD_LVLT   interface dpeths on love bielak model, left, (km)
!     HDD_LVRT   interface dpeths on love bielak model, right, (km)
!     HDD_RLLT   interface dpeths on rayleigh bielak model, left, (km)
!     HDD_RLRT   interface dpeths on rayleigh bielak model, right, (km)
!     RH1D_LT    densities on left Bielak side (kg/m**3)
!     RH1D_RT    densities on right Bielak side (kg/m**3)
!     ROD_LVLT   unflattened love density left density (g/cm**3) 
!     ROD_LVRT   unflattened love density right density (g/cm**3)
!     ROD_RLLT   unflattened rayleigh density left density (g/cm**3) 
!     ROD_RLRT   unflattened rayleigh density right density (g/cm**3)
!     VP1D_LT    Vp 1D model on left Bielak side (m/s)
!     VP1D_RT    Vp 1D model on right Bielak side (m/s)
!     VPD_LVLT   unflattened love vp left model (km/s)
!     VPD_LVRT   unflattened love vp right model (km/s)
!     VPD_RLLT   unflattened rayleigh vp left model (km/s)
!     VPD_RLRT   unflattened rayleigh vp right model (km/s)
!     VS1D_LT    Vs 1D model on left Bielak side (m/s)
!     VS1D_RT    Vs 1D model on right Biealk side (m/s)
!     VSD_LVLT   unflattened love vs left model (km/s)
!     VSD_LVRT   unflattened love vs right model (km/s)
!     VSD_RLLT   unflattened rayleigh vs left model (km/s)
!     VSD_RLRT   unflattened rayleigh vs right model (km/s)
!     Z1D_LT     interface depths on Bielak left side (m)
!     Z1D_RT     interface depths on Bielak right side (m)
!
!.... local variables
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH 
      TYPE (MOD1D_INFO) MOD1D
!.... local variables
      REAL*8, ALLOCATABLE :: VP1D(:), VS1D(:), RH1D(:), Z1D(:), QP1D(:), QS1D(:)
      INTEGER*4 NWORK,NL1D,ICNIBLK
      LOGICAL*4 LFLIP, LVERB
      PARAMETER(LFLIP = .TRUE.) !flip model
      PARAMETER(LVERB = .TRUE.) !can output unflattened 1D models  
!
!----------------------------------------------------------------------------------------!
!
!.... fill left model
      NWORK = (MSH%NLXI - 1)*ICNIBLK(MSH%NNPG,MSH%CNNPG) + MSH%NLXI
      ALLOCATE(VP1D(NWORK))
      ALLOCATE(VS1D(NWORK))
      ALLOCATE(RH1D(NWORK))
      ALLOCATE(Z1D (NWORK))
      ALLOCATE(QP1D(NWORK)) 
      ALLOCATE(QS1D(NWORK)) 
      CALL G1DBLK(msh%NNPG,NGNOD, msh%NNPG,NWORK,msh%NELEM,msh%NLXI, NGNOD, 'L',    &
                  msh%CDOMAIN,msh%CNNPG, msh%IENG, &
                  msh%XIPTS,msh%XLOCS,msh%ZLOCS,msh%DENS,   &
                  msh%QP,msh%QS,msh%ECOEFF,    &
                  NL1D,msh%XBLKL,msh%XMOD0, Z1D,VP1D,VS1D,RH1D,QP1D,QS1D)
      mod1d%NL1D_LT = NL1D
      ALLOCATE(mod1d%VP1D_LT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%VS1D_LT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%RH1D_LT(mod1d%NL1D_LT))
      ALLOCATE(mod1d% Z1D_LT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%QP1D_LT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%QS1D_LT(mod1d%NL1D_LT))
      mod1d%VP1D_LT(1:MOD1D%NL1D_LT) = VP1D(1:NL1D)
      mod1d%VS1D_LT(1:MOD1D%NL1D_LT) = VS1D(1:NL1D)
      mod1d%RH1D_LT(1:MOD1D%NL1D_LT) = RH1D(1:NL1D)
      mod1d% Z1D_LT(1:MOD1D%NL1D_LT)  = Z1D(1:NL1D)
      mod1d%QP1D_LT(1:mod1d%NL1D_LT) = QP1D(1:NL1D)
      mod1d%QS1D_LT(1:mod1d%NL1D_LT) = QS1D(1:NL1D)
!
!.... rayleigh and love wave models
      ALLOCATE(mod1d%VPD_RLLT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%VSD_RLLT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%ROD_RLLT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%HDD_RLLT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%VPD_LVLT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%VSD_LVLT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%ROD_LVLT(mod1d%NL1D_LT))
      ALLOCATE(mod1d%HDD_LVLT(mod1d%NL1D_LT))
      CALL MKERTHMODS(MOD1D%NL1D_LT,LFLIP,LVERB, Z1D,VP1D,VS1D,RH1D,   &   
                      MOD1D%VPD_RLLT,MOD1D%VSD_RLLT,MOD1D%ROD_RLLT,MOD1D%HDD_RLLT,  &
                      MOD1D%VPD_LVLT,MOD1D%VSD_LVLT,MOD1D%ROD_LVLT,MOD1D%HDD_LVLT)
!
!.... fill right model
      CALL G1DBLK(msh%NNPG,NGNOD, msh%NNPG,NWORK,msh%NELEM,msh%NLXI, NGNOD, 'R',    &
                  msh%CDOMAIN,msh%CNNPG, msh%IENG, &
                  msh%XIPTS,msh%XLOCS,msh%ZLOCS,msh%DENS,  &
                  msh%QP,msh%QS,msh%ECOEFF,    &
                  NL1D,msh%XBLKR,msh%XMOD1, Z1D,VP1D,VS1D,RH1D,QP1D,QS1D)
      mod1d%NL1D_RT = NL1D
      ALLOCATE(mod1d%VP1D_RT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%VS1D_RT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%RH1D_RT(mod1d%NL1D_RT))
      ALLOCATE(mod1D% Z1D_RT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%QP1D_RT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%QS1D_RT(mod1d%NL1D_RT))
      MOD1D%VP1D_RT(1:MOD1D%NL1D_RT) = VP1D(1:NL1D)
      MOD1D%VS1D_RT(1:MOD1D%NL1D_RT) = VS1D(1:NL1D)
      MOD1D%RH1D_RT(1:MOD1D%NL1D_RT) = RH1D(1:NL1D)
      MOD1D% Z1D_RT(1:MOD1D%NL1D_RT)  = Z1D(1:NL1D)
      mod1d%QP1D_RT(1:mod1d%NL1D_RT) = QP1D(1:NL1D)
      mod1d%QS1D_RT(1:mod1d%NL1D_RT) = QS1D(1:NL1D)
!
!.... rayleigh and love wave models
      ALLOCATE(mod1d%VPD_RLRT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%VSD_RLRT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%ROD_RLRT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%HDD_RLRT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%VPD_LVRT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%VSD_LVRT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%ROD_LVRT(mod1d%NL1D_RT))
      ALLOCATE(mod1d%HDD_LVRT(mod1d%NL1D_RT))
      CALL MKERTHMODS(MOD1D%NL1D_RT,LFLIP,LVERB, Z1D,VP1D,VS1D,RH1D,   &   
                      MOD1D%VPD_RLRT,MOD1D%VSD_RLRT,MOD1D%ROD_RLRT,MOD1D%HDD_RLRT,  &
                      MOD1D%VPD_LVRT,MOD1D%VSD_LVRT,MOD1D%ROD_LVRT,MOD1D%HDD_LVRT) 
!
!.... clean space
      DEALLOCATE(VP1D)
      DEALLOCATE(VS1D)
      DEALLOCATE(RH1D)
      DEALLOCATE(Z1D)
      DEALLOCATE(QP1D)
      DEALLOCATE(QS1D) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FILL_SRCPRM(IFREQ, VFAST,AZTOL,AOITOL, SRC) 
!
!     Fills the source permutation list for this frequency group
!
!     INPUT      MEANING
!     -----      ------- 
!     AOITOL     angle of incidence tolerance (degreess)
!     AZTOL      azimuth tolerance (degrees)
!     IFREQ      frequency we are modeling
!     VFAST      max velocity at base of models
!  
!     OUTPUT     MEANING
!     ------     -------
!     ISGPTR     source group pointer
!     ISRCPRM    source permuatation list
!     NSG        number of source groups
!     PYAVG      average py for source group
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (SRC_INFO) SRC 
      REAL*8, INTENT(IN) :: VFAST, AZTOL, AOITOL 
      INTEGER*4, INTENT(IN) :: IFREQ
      REAL*8, ALLOCATABLE :: PYAVG(:) 
      INTEGER*4, ALLOCATABLE :: ISGPTR(:), ISRCPRM(:)  
      LOGICAL*4 LVERB
!
!----------------------------------------------------------------------------------------!
!
!.... set space 
      ALLOCATE(PYAVG(SRC%NSRC))
      ALLOCATE(ISGPTR(SRC%NSRC+1))
      ALLOCATE(ISRCPRM(SRC%NSRC))
      PYAVG(1:SRC%NSRC) = 0.D0
      ISGPTR(1:SRC%NSRC) = 0 
      ISRCPRM(1:SRC%NSRC) = 0 
      LVERB = .FALSE. 
      IF (IFREQ == 1) LVERB = .TRUE.
!.... generate the source list
      CALL SRCLIST4(SRC%NSRC,LVERB, VFAST,AZTOL,AOITOL, SRC%PYTAB(IFREQ,1:SRC%NSRC), &
                    SRC%NSG,ISGPTR,ISRCPRM,PYAVG) 
!.... resize 
      ALLOCATE(SRC%ISGPTR(SRC%NSG+1))
      SRC%ISGPTR(1:SRC%NSG+1) = ISGPTR(1:SRC%NSG+1)
      ALLOCATE(SRC%ISRCPRM(SRC%NSRC))
      SRC%ISRCPRM(1:SRC%NSRC) = ISRCPRM(1:SRC%NSRC)
      ALLOCATE(SRC%PYAVG(SRC%NSG))
      SRC%PYAVG(1:SRC%NSG) = PYAVG(1:SRC%NSG)
!
!.... clean
      DEALLOCATE(ISGPTR)
      DEALLOCATE(ISRCPRM)
      DEALLOCATE(PYAVG) 
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FILL_SRCPRM_SB(IFREQ,LSURF, VFAST,AZTOL,AOITOL, SRC)
!
!     Fills the source permutation list for this frequency group for 
!     surface wave modeling or body wave modeling
!
!     INPUT      MEANING
!     -----      ------- 
!     AOITOL     angle of incidence tolerance (degreess)
!     AZTOL      azimuth tolerance (degrees)
!     IFREQ      frequency we are modeling
!     LSURF      True -> want surface wave list, False -> body wave
!     VFAST      max velocity at base of models
!  
!     OUTPUT     MEANING
!     ------     -------
!     ISGPTR     source group pointer
!     ISRCPRM    source permuatation list
!     NSG        number of source groups
!     PYAVG      average py for source group
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (SRC_INFO) SRC 
      REAL*8, INTENT(IN) :: VFAST, AZTOL, AOITOL 
      INTEGER*4, INTENT(IN) :: IFREQ
      LOGICAL*4, INTENT(IN) :: LSURF
      REAL*8, ALLOCATABLE :: PYAVG(:), PYWORK(:) 
      INTEGER*4, ALLOCATABLE :: ISGPTR(:), ISRCPRM(:), IMASK(:)  
      INTEGER*4 NSRC, ISRC, JSRC
      LOGICAL*4 LVERB, LBODY 
!
!----------------------------------------------------------------------------------------!
!
!.... intial check
      IF (LSURF) THEN
         NSRC = src%NSRC_SRF
         LBODY = .FALSE.
      ELSE
         NSRC = src%NSRC_BDY
         LBODY = .TRUE.
      ENDIF
      IF (NSRC == 0) THEN
         src%NSG = 0
         RETURN
      ENDIF
      ALLOCATE(PYWORK(NSRC)) 
      ALLOCATE(IMASK(NSRC)) 
      JSRC = 0
      DO 1 ISRC=1,src%NSRC
         IF (LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
             LBODY .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
            JSRC = JSRC + 1
            PYWORK(JSRC) = src%PYTAB(IFREQ,ISRC)
            IMASK(JSRC) = ISRC
         ENDIF
    1 CONTINUE  
      IF (JSRC /= NSRC) THEN
         WRITE(*,*) 'fill_srcprm_sb: Warning jsrc /= nsrc!'
      ENDIF 
!
!.... generate the source list
      ALLOCATE(ISGPTR(NSRC+1))
      ALLOCATE(ISRCPRM(NSRC))
      ALLOCATE(PYAVG(NSRC)) 
      LVERB = .FALSE.
      IF (IFREQ == 1) LVERB = .TRUE.
      CALL SRCLIST4(NSRC,LVERB, VFAST,AZTOL,AOITOL, PYWORK, & 
                    src%NSG,ISGPTR,ISRCPRM,PYAVG) 
!.... resize 
      ALLOCATE(src%ISGPTR(src%NSG+1))
      src%ISGPTR(1:src%NSG+1) = ISGPTR(1:src%NSG+1)
      !CALL ICOPY(src%NSG+1,ISGPTR,1,src%ISGPTR,1)
      ALLOCATE(src%ISRCPRM(NSRC))
      DO 2 ISRC=1,NSRC
         src%ISRCPRM(ISRC) = IMASK(ISRCPRM(ISRC))
    2 CONTINUE
      ALLOCATE(src%PYAVG(src%NSG))
      CALL DCOPY(src%NSG,PYAVG,1,src%PYAVG,1)
!
!.... clean
      DEALLOCATE(ISGPTR)
      DEALLOCATE(ISRCPRM)
      DEALLOCATE(PYAVG)
      DEALLOCATE(PYWORK)
      DEALLOCATE(IMASK) 
      RETURN
      END


