      subroutine mkerthmods(nl,lflip,lverb, z1d,vp1d,vs1d,rh1d, &
                            vpd_rl,vsd_rl,rod_rl,hdd_rl, vpd_lv,vsd_lv,rod_lv,hdd_lv)
!
!
!     Unflattens the 1D velocity model for surface wave calculation
!
!     input      meaning
!     -----      ------- 
!     lflip      True -> flip coordinate system so z is + up  
!     lverb      controls verbosity level
!     nl         number of interfaces
!     rh1d       1D density model (kg/m**3)
!     vp1d       1D compressional velocity model (m/s)
!     vs1d       1D shear velocity model (m/s)
!     z1d        1D depths of interfaces z locations
!     
!     output     meaning
!     ------     ------- 
!     h1d_lv     Unflattened interface depths in Love model (km)
!     h1d_rl     Unflattened interface depths in Rayleigh model (km)
!     rod_lv     Unflattened Love density model (g/cm**3) 
!     vpd_lv     Unflattened Love p velocity model (km/s)
!     vsd_lv     Unflattened Love s velocity model (km/s)
!     rod_rl     Unflattened Rayleigh density model (g/cm**3) 
!     vpd_rl     Unflattened Rayleigh p velocity model (km/s)
!     vsd_rl     Unflattened Rayleigh s velocity model (km/s)
!
!.... variable declarations
      real*8, intent(in) :: z1d(nl), vp1d(nl), vs1d(nl), rh1d(nl) 
      logical*4, intent(in) :: lflip, lverb
      real*8, intent(out) :: vpd_rl(nl), vsd_rl(nl), rod_rl(nl), hdd_rl(nl),  &
                             vpd_lv(nl), vsd_lv(nl), rod_lv(nl), hdd_lv(nl)  
!.... local variables
      real*8 ztop
      real*8, allocatable :: h1d(:)
!
!----------------------------------------------------------------------------------------!
!
!.... set space and convert 
      allocate(h1d(nl+1))
      h1d(1:nl+1) = 0.d0
      ztop = z1d(1) 
      do i=1,nl 
         if (lflip) then
            h1d(i) = ztop - z1d(i) 
         else
            h1d(i) = z1d(i)
         endif
         vpd_rl(i) = vp1d(i)*0.001d0 !m/s -> km/s
         vsd_rl(i) = vs1d(i)*0.001d0 !m/s -> km/s
         rod_rl(i) = rh1d(i)*0.001d0 !kg/m^3 -> g/cm^3
         vpd_lv(i) = vpd_rl(i) !love vp 
         vsd_lv(i) = vsd_rl(i) !love vs 
         rod_lv(i) = rod_rl(i) !love density
      enddo
      h1d(nl+1) = h1d(nl)
      do i=1,nl
         hdd_rl(i) = (h1d(i+1) - h1d(i))*0.001d0 !m -> kg
         hdd_lv(i) = hdd_rl(i) !layers thicknesses
      enddo
!
!.... unflatten 
!     call unflat(1, vpd_rl,vsd_rl, rod_rl, hdd_rl, nl) !rayleigh waves
!     call unflat(2, vpd_lv,vsd_lv, rod_lv, hdd_lv, nl) !love 

      if (lverb) then
         write(*,*) ''
         write(*,*) 'mkerthmods: 1D Input Model Summary'
         write(*,905)
         do i=1,nl
            write(*,907) i,vp1d(i),vs1d(i),rh1d(i),h1d(i)*1.d-3
         enddo
         write(*,*) 
         write(*,*) 'mkerthmods: Unflattened Rayleigh Wave Summary'
         write(*,906)
         do i=1,nl
            write(*,907) i, vpd_rl(i), vsd_rl(i), rod_rl(i), hdd_rl(i)
         enddo
         write(*,*)
         write(*,*) 'mkerthmods: Unflattened Love Wave Summary'
         write(*,906)
         do i=1,nl
            write(*,907) i, vpd_lv(i), vsd_lv(i), rod_lv(i), hdd_lv(i) 
         enddo
      endif
      deallocate(h1d)
  905 format('             Layer   Vp (m/s)   Vs (m/s)  Dens(kg/m^3)       Depth (km)',/,&
             '             -----   --------   --------  ------------       ----------')
  906 format('             Layer   Vp(km/s)   Vs(km/s)  Dens(g/cm^3)   Thickness (km)',/,&
             '             -----   --------   --------  ------------   --------------')
  907 format(12x,i6,2x,f9.4,2x,f9.4,3x,f9.4,7x,f12.4) 
      if (lverb) write(*,*) 
      return
      end 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      subroutine unflat(jcom, vp, vs, ro, h, n)
!
!     This subroutine undoes the earth-flattening corrections
!     This adapted from earthsubs.f. For reference see biswas 1972.
!
!     jcom    = 1 rayleigh waves ; <> 1 love waves
!     h(i)    = thickness of layer i
!     ro(i)   = density of layer i
!     vp(i)   = P wavespeed in layer i
!     vs(i)   = S wavespeed in layer i
!
!     h, ro, vp, and vs are replace by their flattened vaules on return
!
      implicit real*8(a-h,o-z)
      real*8, intent(inout) :: vp(n), vs(n), ro(n), h(n)
      integer*4, intent(in) :: jcom, n
      parameter(a = 6371.0008d0)
!
!----------------------------------------------------------------------------------------!
!
      pwr = 2.275d0
      if (abs(jcom) > 1) pwr = 5.d0
      atp = a**pwr
!
!.... unflatten model
      nm = n - 1 
      hs = 0.d0 
      do 1 i = 1, n
         ht   = hs
         hs   = hs + h(i)
         h(i) = ht
    1 continue

      h(1) = a
      do 2 i = 2,n 
        h(i) = a*dexp(-h(i)/a)
    2 continue

      do 3 i = 1, nm
         ii = i + 1 
         fltd  = dlog(h(i)/h(ii))
         dif   = (1.d0/h(ii) - 1.d0/h(i))*a/fltd
         difr  = h(i)**pwr - h(ii)**pwr
         ro(i) = ro(i)*(fltd*(a**pwr)*pwr)/difr
         vp(i) = vp(i)/dif
         vs(i) = vs(i)/dif
    3 continue
!.... half space scaling
      fact = a/h(n)
      facti = (h(n)**pwr)/atp
      vp(n) = vp(n)/fact
      vs(n) = vs(n)/fact
      ro(n) = ro(n)/facti
      do 4 i = 1, nm
         h(i) = h(i) - h(i+1)
    4 continue 
      h(n) = 0.d0
      return
      end

