        subroutine mkhomogsr(idfl,idfl2,c1,s1,c2,s2,modeuse, lverb,lbin)

c  This subroutine reads in output from the excitation program (earth.f).
c  Equations are taken from mendiguren (1977) but the signs of
c  imaginary parts are opposite to those given.  The radial component
c  also is multiplied by -1.  greens functions have been tested
c  against those computed from a normal mode sum for azimuths in all
c  four quadrants.  differences with mendiguren possibly stem from the
c  different coordinate system used in the computation of the excitations
c  and ambiguities in mendiguren's paper.  This program does agree with Aki
c  and Richards derivation (1981).

        implicit real*8 (a-h,o-z)
        implicit integer*4(i-n)
        include 'sizes.inc'
        include 'commons.inc'
        !include 'units.inc'
       integer*4 iinf1,iinf2,iinf3,iinf4,iinf5
           parameter (iinf1=15,iinf2=16,iinf3=19,iinf4=21,iinf5=23)
       integer*4 iouf1,iouf2,iouf3,iouf4
           parameter (iouf1=17,iouf2=18,iouf3=20,iouf4=22)

c       real*8 sddep,w,ccc,gmmm,py1,py2,py3,ddt
        real*8 sddep,w,ccc,gmmm,py1,py2,py3,ddf
        real*8 c1,c2,s1,s2,b1,b2,gb1,gb2
        logical*4 lverb, lbin
          
        common/strain/strn(lgrm,9)

        complex*16 u,wu,dudr,dudz,dwdr,dwdz,erz,gr
        complex*16 v,dvdr,dvdz,erp,gl
        dimension sddep(lsd),py(3)

        data pi/3.1415926535897931d0/
        parameter(eps = 1.11d-7)

        c2p1 = 0.5d0*(1.d0+c2)
        c2m1 = 0.5d0*(1.d0-c2)
        s22 = 0.5d0*s2
           
c  open files containing excitations
c       read(idfl) nsrce,npts,ddt,jcom,nbran
        if (lbin) then
           read(idfl) nsrce,nom,ddf,jcom,nbran
        else
           read(idfl,*) nsrce,nom,ddf,jcom,nbran
        endif

c-- modeuse is a counter to let us know if any modes are actually summed here
        modeuse = 0
          
c  check that sample intervals are consistent
        !if (df.ne.sngl(ddf)) then
        if (dabs(df - ddf).gt.eps) then
          df = sngl(ddf)
          write(*,103) df
103       format(/,' Frequency interval changed to ',f5.2)
        endif

c  check that component and input files are consistent
        if((jcom.ne.1)) then
          write(*,'(a)') ' ***WRONG FILE TYPE***'
          nom = 0
          close(idfl)
          return
        endif
          
c  choose source depth
        if (lbin) then
           read(idfl) (sddep(i),i = 1,nsrce)
        else
           read(idfl,*) (sddep(i),i = 1,nsrce)
        endif

        do  i = 1,nsrce
          sdep(i) = sddep(i)
        enddo
c  find nearest source depth in file to event depth
        !call srcget(d0,nsdp)
        call srcget
        sourcd = sdep(nsdp)
c initialize arrays and constants

c       nper  = npts/2
        nper  = nom
        nh1   = nper+1
        npts = 2*nom
        npts2 = npts+2
c       trec  = npts*dt
        trec  = 1.d0/df
c       df    = 2.*pi/trec
        dom   = 2.d0*pi*df
c       nlen  = npts
        nlen  = nom
 
        do j = 1,9
          do i = 1,npts2
            strn(i,j) = 0.d0
          enddo
        enddo
          
c calculate instrument response: this fills up the resp.  We can make this an
c optional read at the main program, but for now we just set the response to "1"
c       call inresp(npts)
        do i = 2, nh1
          resp(i) = dcmplx(1.d0,0.d0)
        enddo

c---disply asks for input on number of modes and source depth; skip this now
c       call disply

c   fact accounts for distance, sign conventions, and fft normalization
c   fft sign conventions are as in mendiguren
c       fact = -fmom*sqrt(1.e3/dist)*2./trec
c---it appears that the fft used in mkhomo scales by 1/npts.  The actual scaling
c   shoud be df/nper, and as nper = npts/2, then we have 2*df/npts or (2/trec)/npts.
c   So, since it was already scaled by npts, we are left with 2/trec.   We eliminate
c   this below as we do not fft.
c
        fact = -fmom*dsqrt(1.d3/dist)

c  read in vertical and radial excitation functions
        do 8888 kk = 1,nbran
          if (lbin) then
             read(idfl,end = 10) nb
          else
             read(idfl,*,end = 10) nb
          endif
          do 888 j = 1,nh1
c read in rayleigh excitations
            iki = 2*(nper+2-j)
            ikr = iki-1
            if (lbin) then
               read(idfl,err = 222) w,ccc,gmmm,b1,b2,gb1,gb2
            else
               read(idfl,*,end=222) w,ccc,gmmm,b1,b2,gb1,gb2
            endif
            c = ccc !sngl(ccc)
            gamm = gmmm !sngl(gmmm)
            do 43 k = 1,nsrce
              if (lbin) then
                 read(idfl) py1,py2,py3
              else
                 read(idfl,*) py1,py2,py3
              endif
              if (k.eq.nsdp) then
                py(1) = py1 !sngl(py1)
                py(2) = py2 !sngl(py2)
                py(3) = py3 !sngl(py3)
              endif
43          continue
            if (nb.ge.mdmin.and.nb.le.mdmx) then
              modeuse = modeuse + 1
              xk = w/c
              wvd = xk*dist
              cwd = dcos(wvd)
              swd = dsin(wvd)
                
              gr = dcmplx( py(1)*(s2*sol(6)+c2m1*sol(3)+
     &             c2p1*sol(2))+py(3)*sol(1),py(2)
     ;            *(c1*sol(4)+s1*sol(5)))
         
              atn = dexp(-gamm*dist)
        
              u  = gr*dcmplx(atn*b2,0.d0)*dcmplx(swd, cwd)*1.d-5
              wu = gr*dcmplx(atn*b1,0.d0)*dcmplx(cwd,-swd)*1.d-5
                
              dudr = -dcmplx(0.d0,xk)*u
              dudz = -gb2*u/b2
                
              dwdr = -dcmplx(0.d0,xk)*wu
              dwdz = -gb1*wu/b1
                
              erz = .5d0*(dwdr+dudz)

c------
c       strn(*,7)       w (vertical)  component
c       strn(*,8)       u (radial)    component
c       strn(*,9)       v (transverse)component
c------
        
              if (ivel.gt.0) then
                strn(ikr,7) = strn(ikr,7)+ dimag(wu*w)
                strn(iki,7) = strn(iki,7)- real(wu*w)
                strn(ikr,8) = strn(ikr,8)+ dimag(u*w)
                strn(iki,8) = strn(iki,8)- real(u*w)
              else
                strn(ikr,7) = strn(ikr,7)+ real(wu)*1.d5
                strn(iki,7) = strn(iki,7)+ dimag(wu)*1.d5
                strn(ikr,8) = strn(ikr,8)+ real(u)*1.d5
                strn(iki,8) = strn(iki,8)+ dimag(u)*1.d5
              endif
                
              strn(ikr,1) = strn(ikr,1)+ real(dudr)
              strn(iki,1) = strn(iki,1)+ imag(dudr)

              strn(ikr,3) = strn(ikr,3)+ real(dwdz)
              strn(iki,3) = strn(iki,3)+ imag(dwdz)

              strn(ikr,5) = strn(ikr,5)+ real(erz)
              strn(iki,5) = strn(iki,5)+ imag(erz)

            endif
888       continue
222       if (lverb) write(*,333) nb
          if (nb.eq.mdmx) go to 11

8888    continue
c11     close(idfl)
c
!11      idfl = idfl+2
 11     continue
        if (lverb) write(*,*) ' Finished Rayleigh Waves'
c  read in love excitations
c       read(idfl) nsrce,npts,ddt,jcom,nbran
        if (lbin) then
           read(idfl2) nsrce,nom,ddf,jcom,nbran
        else
           read(idfl2,*) nsrce,nom,ddf,jcom,nbran
        endif

c  check that sample intervals and npts are consistent
        !if (df.ne.sngl(ddf)) then
        if (dabs(df - ddf).gt.eps) then
          write(*,'(a)') ' Error: Love file df is inconsistent!'
          write(*,*) ' df, ddf = ', df, dff
          npts = 0
          close(idfl2)
          return
        endif

        if (nom.ne.nlen) then
          write(*,'(a)') ' Error: Love file nom is inconsistent!'
          nom = 0
          close(idfl2)
          return
        endif
c  check that component and input files are consistent
        if((jcom.ne.2)) then
          write(*,'(a)') ' ***WRONG FILE TYPE***'
          nom = 0
          close(idfl2)
          return
        endif

c  choose source depth
        if (lbin) then
           read(idfl2) (sddep(i),i = 1,nsrce)
        else
           read(idfl2,*) (sddep(i),i = 1,nsrce)
        endif
        do  i = 1,nsrce
          sdep(i) = sddep(i)
        enddo

c  find nearest source depth in file to event depth
        call srcget
        !if(sourcd.ne.sdep(nsdp)) then
        if (dabs(sourcd - sdep(nsdp)).gt.eps) then
          write(*,'(a)') ' Error: Love file source depths inconsistent!'
          npts = 0
          close(idfl2)
          return
        endif

c---loop over all available modes
        do 7777 kk = 1,nbran
          if (lbin) then
             read(idfl2,end = 10) nb
          else
             read(idfl2,*,end = 10) nb
          endif
          do 777 j = 1,nh1
            iki = 2*(nper+2-j)
            ikr = iki-1
c  read in love excitations
            if (lbin) then
               read(idfl2,err = 555) w,ccc,gmmm,b1,gb1
            else
               read(idfl2,*,end = 555) w,ccc,gmmm,b1,gb1
            endif
            c = ccc
            gamm = gmmm
            do k = 1,nsrce
              if (lbin) then
                 read(idfl2) py1,py2
              else
                 read(idfl2,*) py1,py2
              endif
              if (k.eq.nsdp) then
                py(1) = py1
                py(2) = py2
              endif
            enddo

            if (nb.ge.mdmin.and.nb.le.mdmx) then

              xk = w/c
              wvd = xk*dist
              cwd = dcos(wvd)
              swd = dsin(wvd)
                
              gl = dcmplx( py(2)*(-c1*sol(5)+s1*sol(4)),
     &           py(1)*(c2*sol(6)+s22*(sol(3)-sol(2))) )
                        
              atn = dexp(-gamm*dist)
                                        
              v = gl*dcmplx(atn*b1,0.d0)*dcmplx(cwd,-swd)*1.d-5
                
              dvdr = -dcmplx(0.d0,xk)*v
              dvdz = 0.d0
              if(b1.ne.0.d0) dvdz = -.5d0*gb1*v/b1

              if(ivel.gt.0) then
                strn(ikr,9) = strn(ikr,9)+ dimag(v*w)
                strn(iki,9) = strn(iki,9)- dreal(v*w)
              else
                strn(ikr,9) = strn(ikr,9)+ dreal(v)*1.d5
                strn(iki,9) = strn(iki,9)+ dimag(v)*1.d5
              endif
                                                
              strn(ikr,2) = 0.d0
              strn(iki,2) = 0.d0

              strn(ikr,4) = strn(ikr,4)+ dreal(dvdz)
              strn(iki,4) = strn(iki,4)+ dimag(dvdz)

              erp = .5d0*dvdr
              strn(ikr,6) = strn(ikr,6)+ dreal(erp)
              strn(iki,6) = strn(iki,6)+ dimag(erp)

            endif

777       continue
555       if (lverb) write(*,333) nb
          if(nb.eq.mdmx) go to 10

7777    continue

c10     close(idfl)
10      if (lverb) write(*,*) ' Finished Love Waves'
333     format(' mode ',i2,' completed')  
  
c   include travel time,pi/4 phase shifts and instr. response

        do 666 kk = 2,nh1
           iki = 2*kk
           ikr = iki-1
           w = dfloat(kk-1)*dom
           pht = pi*.25d0+w*to
           cph = dcos(pht)
           sph = dsin(pht)
           respr = dreal(resp(kk))
           respi = dimag(resp(kk))
           crl = (respr*cph-respi*sph)*fact
           cri = (respr*sph+respi*cph)*fact
           do jj = 1,9
             save = (strn(ikr,jj)*crl-strn(iki,jj)*cri)
             strn(iki,jj) = (strn(ikr,jj)*cri+strn(iki,jj)*crl)
             strn(ikr,jj) = save
           enddo
666     continue
  
c fft to the time domain
c       do j = 1,9
c         call fftl(strn(1,j),nlen,-2,ier)
c         if(ier.ne.0) then
c           write(*,'(a)') ' fft error!'
c           npts = 0
c           return
c         endif
c       enddo  

        return
        end

      subroutine srcget
c subroutine finds source depth in input file nearest to desired depth
c       include 
c     +  '/home/keith/codes/surface_waves/earth/common_files/sizes.inc'
c       include 
c     +  '/home/keith/codes/surface_waves/earth/common_files/commons.inc'
     
       include 'sizes.inc'
       include 'commons.inc'
      real*8 dnear, dif

      dnear=9999.d0
      do 9 i=1,nsrce
      dif=dabs(d0-sdep(i))
      if(dif.gt.dnear) go to 9
      nsdp=i
      dnear=dif
    9 continue
      return
      end 

