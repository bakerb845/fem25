c----------------------------------------------------------------------
c
c       Compute the Green's functions for a 1D medium using TH (Thomson-Haskell).  
c
c       In this version we use the model contained in vp1disc and vs1disc
c       which nonmialy is the 1D model corrected for dipsersion (i.e., they
c       are the wavespeeds in the discrete model, as opposed to the continuous
c       model.
c
c       This version differs from the orignal haskgrn only in the use of
c       vp1disc and vs1disc instead of vp1d and vs1d.
c
c       This routine is the computative "engine" of the T-H method.
c       It impliments the concepts in Haskell (1953; H53) and Haskell
c       (1962; H62) to compute the Green's functions in layered media
c       do to an external plane wave source of unit amplitude
c
c       Arguments:
c       df      The frequency 
c       c       The apparent velocity
c       ro      The density array
c       h       The thickness array
c
c       This version for 1 frequency = df
c
c       Author:  S. Roecker 6/09
c
        subroutine haskgrnsh(nl,nzpts,lflip, df,ain, z1d,
     ;                       vs1d,rh1d,depth, ugrn1f,wgrn1f, ierr)  
!
!       input      meaning
!       -----      ------- 
!       ain        angle of incidence (degrees) 
!       depth      depths at which to evaluated 1d greens fns (m) 
!       df         frequency (hz) 
!       lflip      true -> flip then flip coordinate system 
!       nl         number of interfaces (nlayers + 1)  
!       nzpts      number of points in z vectors
!       rh1d       density at interfaces (kg/m**3) 
!       vs1d       S velocity at interfaces (m/s) 
!       z1d        interfaces depths (m)  
! 
!       output     meaning
!       ------     ------- 
!       ugrn1f     1d greens functions in horizontal  
!       wgrn1f     1d greens functions in vertical (+z down) 

        implicit none

        real*8, intent(in) :: z1d(nl),vs1d(nl),rh1d(nl),
     ;                        depth(nzpts),ain,df
        integer*4, intent(in) :: nl,nzpts
        logical*4, intent(in) :: lflip
        complex*16, intent(out) :: ugrn1f(nzpts),wgrn1f(nzpts)
        integer*4, intent(out) :: ierr
c       include 'dimension.inc'
c       include 'common.inc'
c       include 'init.inc'

        real*8, allocatable :: r2(:), r2m(:), e2(:), 
     ;                         s02(:), c02(:), h1d(:), h(:) 

        complex*16 quresp
        real*8  x,ur,ui,pic,hc, b11,b21, d11,d21, a11,a12,a21,a22,
     ;          ztop, covb, dept, dloc, qiu, qru
        real*8 z02
        real*8 cv
        real*8 mum

        integer*4 i, m, iz, izb,ize,idir, ibov, ibovo, ifin

        real*8 pi180,pi
        parameter(pi180 = 0.017453292519943295d0,
     ;            pi = 3.1415926535897932384626434d0)

!...... horizontal wavespeed
        IERR = 0
        cv = vs1d(nl)/dsin(ain*pi180)

        pic = 2.d0*pi*df/cv
c       pic = df/cv

c       write(*,*) ' cv, pic = ', cv, pic

!
!...... may need to flip coordinate system in z
        allocate(h(nl-1))
        h(1:nl-1) = 0.d0
        allocate(h1d(nl))
        h1d(1:nl) = 0.d0
        ztop = z1d(1)
        do i=1,nl
           if (lflip) then
              h1d(i) = ztop - z1d(i) 
           else
              h1d(i) = z1d(i)
           endif
        enddo
!
!...... thickness array where layer 1 is the model basement 
        do i = 1, nl-1
          h(i) = h1d(i+1) - h1d(i)
        enddo
!
!...... set space
        allocate(r2(nl+1))
        allocate(r2m(nl+1))
        allocate(e2(nl+1))
        allocate(c02(nl+1))
        allocate(s02(nl+1))
!
!...... Loop over the model, setting up layer specific values
        do m = 1,nl

c--- c/beta
          covb = cv/vs1d(m)
c----r2 is  r(beta)
          x = covb*covb - 1.d0
          r2(m) = dsign(dsqrt(dabs(x)),x)

          if(m.lt.nl) then
            hc = pic*h(m)

c----z02 is Qm in H53
c----c02 and s02 are cos(Qm) and sin(Qm)
            z02 = hc*r2(m)
c----Case of r(beta) being imaginary
            if(cv.gt.vs1d(m)) then
              e2(m) = 1.d0
              c02(m) = dcos(z02)
              s02(m) = dsin(z02)
            else
c----for imaginary case, use hyperbolic sin and cosine
              x = dexp(z02)
              e2(m) = -1.d0
              c02(m) = 0.5d0*(x + 1.d0/x)
              s02(m) = 0.5d0*(x - 1.d0/x)
            end if
          endif

c---define r2m to be mu*r(beta)
          mum = vs1d(m)*vs1d(m)*rh1d(m)
          r2m(m) = mum*r2(m)

        enddo

c---Now propagate to the surface to compute the surface displacements

c----The 2 x 1 b matrix will hold the multiplication
c    of all the a matrices.  This initialization allows
c    the first column of the a matrix to be saved
c    the first time around.

        b11 = 1.d0
        b21 = 0.d0

c----Loop over all layers
        do 7 m = 1,nl-1

          a11 = c02(m)
          a12 = s02(m)/r2m(m)
          a21 = e2(m)*s02(m)*r2m(m)
          a22 = c02(m)
c
c       D is the accumulation of the multiplication of the A
c       matrices as we descend through the layers.  Note that
c       because we start at the surface, the stresses are 0 and
c       so the vector we multiply is (u, w, 0, 0).  This means 
c       that for the entire problem we need save only the first
c       two columns of the multiplication.
c
c       The + and - signs below are a consequence of every other
c       term in A being real and imaginary.  Thus, we multiply:
c               |R I||R| 
c               |I R||I|
c       See H53, page 23 for details.
c
          d11 =  a11*b11 - a12*b21
          d21 =  a21*b11 + a22*b21
          b11 = d11
          b21 = d21

7       continue
c

c
c       qru, qiu are the real and imaginary components of horizontal displacement
c
        d11 = b11*r2m(nl)
        d21 = d11*d11 + b21*b21
c---these are the surface displacements
        qru =  2.d0*b11*r2m(nl)*r2m(nl)/d21
        qiu = -2.d0*b21*r2m(nl)/d21

c       write(*,*) ' qru, qiu = ', qru, qiu
c       read(*,*) pause

c---now propagate backwards, looping over the grid points
        ibovo  = 1
        ifin = 0

        if (lflip) then
           izb = nzpts
           ize = 1
           idir =-1
        else
           izb = 1
           ize = nzpts
           idir = 1
        endif

c---copy over response at the surface
        ugrn1f(1:nzpts) = dcmplx(0.d0,0.d0)
        wgrn1f(1:nzpts) = dcmplx(0.d0,0.d0)
        quresp = dcmplx(qru,qiu)
        if (dabs(depth(izb) - z1d(1)).lt.1.11d-7) then
           ugrn1f(izb) = quresp
           wgrn1f(izb) = dcmplx(0.d0,0.d0)  
           if (lflip) izb = izb + idir
        endif

c       do iz = 2, nz
        do 100 iz=izb,ize,idir
          !dept = dz*(iz-1) + zorig
          if (lflip) then
             dept = ztop - depth(iz)
          else
             dept = depth(iz)
          endif
          
c----Loop over the model, setting up layer specific values
c----dep(1) should the top of the model (e.g., zero) so start at dep(2)
          do i = 2,nl
            if (dept.lt.h1d(i)) go to 4
          enddo
          i = nl
4         ibov = i - 1
          dloc = dept - h1d(i-1)
          if (dabs(dloc).lt.1.11d-7) dloc = 0.d0
          if (dloc.lt.0.) then
            write(*,*) 'haskgrnsh: Error in haskgrnsh dloc is < 0! '
            ierr = 1 
            return
          endif
c         write(*,*) ' iz, dept, dloc, ibov = ', iz, dept, dloc, ibov

c----see if we have entered a new layer.  If so, finish off the previous one
          if (ibov.ne.ibovo) then
            ibovo = ibov
            if (ibov.gt.1) then
              ifin = ifin + 1
              if (ifin.gt.nl) then
                 write(*,*) 'haskgrnsh: Dimension error'
                 ierr = 2
                 return
              endif
              hc = pic*h(ifin)
c----z02 is Qm in H53
              z02 = hc*r2(ifin)

c----Case of r(beta) being imaginary
              if (cv.gt.vs1d(m)) then
                e2(ifin) = 1.d0
                c02(ifin) = dcos(z02)
                s02(ifin) = dsin(z02)
              else
c----for imainary case, use hyperbolic sin and cosine
                x = dexp(z02)
                e2(ifin) = -1.d0
                c02(ifin) = 0.5d0*(x + 1.d0/x)
                s02(ifin) = 0.5d0*(x - 1.d0/x)
              end if
            end if
          end if

c----current layer terms
          hc = pic*dloc
c----z02 is  Qm in H53
          z02 = hc*r2(ibov)
c----Case of r(beta) being imaginary
          if(cv.gt.vs1d(ibov)) then
            e2(ibov) = 1.d0
            c02(ibov) = dcos(z02)
            s02(ibov) = dsin(z02)
          else
c----for imaginary case, use hyperbolic sin and cosine
            x = dexp(z02)
            e2(ibov) = -1.d0
            c02(ibov) = 0.5d0*(x + 1.d0/x)
            s02(ibov) = 0.5d0*(x - 1.d0/x)
          end if

c----The 4 x 2 b matrix will hold the multiplication
c    of all the a matrices.  This initialization allows
c    the first two columns of the a matrix to be saved
c    the first time around.

c---this part you only need do when you enter a new layer.  Save
c   results up to ifin then then do one ore iteration for the current layer.

          b11 = 1.d0
          b21 = 0.d0

c----Loop over all layers
          do  m = 1,ibov
            a11 = c02(m)
            a12 = s02(m)/r2m(m)
            a21 = e2(m)*s02(m)*r2m(m)
            a22 = c02(m)
            d11 = a11*b11 - a12*b21
            d21 = a21*b11 + a22*b21
            b11 = d11
            b21 = d21
          enddo

          ur =  b11*qru
          ui =  b11*qiu 

c       write(*,*) ' ur, ui, b11 = ', ur, ui, b11
c       read(*,*) pause

          ugrn1f(iz) = dcmplx(ur, ui)
          wgrn1f(iz) = dcmplx(0.d0, 0.d0)

  100   continue
c---end loop over nodes
!
!...... clean space
        deallocate(r2)
        deallocate(r2m)
        deallocate(e2)
        deallocate(c02)
        deallocate(s02)
        deallocate(h1d)
        deallocate(h) 
        return
        end

