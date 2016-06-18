c------------------------------------------------------------------
c
c   IIR filter program.  Based on the iir.c code provided by PASSCAL.
c       These filters are Butterworth highpass and lowpass filters
c       Both filters may be implemented or only one filter.
c
c       If nh = 0 , then only low pass is implimented
c       If nl = 0 , then only high pass is implimented
c
c       The program to calculate the pole position or to filter does
c       not have any limits.
c
c   Input parameters
c     nh = order of the high pass filter can be 0, or an even number
c               up to 12
c     fh = high pass cutoff frequency
c     nl = order of the low pass filter can be 0, or an even number
c               up to 12
c     fl = low pass cutoff frequency
c     dt = sample rate
c     of = time series to filter
c   Output
c     af = output filtered time series.
c
c------------------------------------------------------------------

        subroutine iir(ofin, afin, dtin, npts, nlin, flin, nhin, fhin)
        implicit real*8 (a-h,o-z)
        real*4 afin(npts), ofin(npts), dtin, flin, fhin
        integer npts, nlin, nhin
c---Local variable definition
        integer i
        complex*16, allocatable :: pl(:), ph(:)
        real*8, allocatable :: af(:), of(:)

c---Nyquist frequency
        dt = dble(dtin)
        fn = 1.d0/(2.d0*dt)

        allocate(af(npts))
        allocate(of(npts))
        of(1:npts) = dble(ofin(1:npts))

c---Temporary output to see if parameters read in correctly
c       write(*,*) ' dt = ', dt
c       write(*,*) ' npts = ', npts
c       write(*,*) ' nl = ', nl
c       write(*,*) ' fl = ', fl
c       write(*,*) ' nh = ', nh
c       write(*,*) ' fh = ', fh
c       write(*,*) ' fn = ', fn


c--check high pass cut-off frequency.  Turn off if cut off is zero, or
c--if cut-off is higher than the Nyquist
        nh = nhin
        fh = dble(fhin)
        if (nh.ne.0) then
          if (fh.eq.0.d0) nh = 0
          if (fh.gt.fn) then
            write(*,*) ' WARNING:  HP freq > Nyquist! '
            nh = 0
          end if
        end if

c---check low pass filter parameters
        nl = nlin
        fl = dble(flin)
        if (nl.ne.0) then
          if(fl.eq.0.d0) then
            nl = 0
          else
            if (fl .gt. fn) then
              write(*,*) ' Warning:  fl greater than Nyquist'
              fl = fn
            end if
          end if
        end if

        if(fl.lt.fh) then
          write(*,*) 'WARNING: LOW PASS CUTOFF IS'
          write(*,*) 'LESS THAN HIGH PASS CUTOFF.'
          write(*,*) 'Low Pass Filter being turned off'
          nl = 0
        end if

c---Temporary output to see if parameters read in correctly
        write(*,*) ' dt = ', dt,dtin
        write(*,*) ' npts = ', npts
        write(*,*) ' nl = ', nl,nlin
        write(*,*) ' fl = ', fl,flin
        write(*,*) ' nh = ', nh,nhin
        write(*,*) ' fh = ', fh,fhin
        write(*,*) ' fn = ', fn

c---Get highpass filter poles if necessary
        b0h = 0.d0
        b0l = 0.d0
        if(nh.ne.0) then
          allocate(ph(nh)) 
          call highpass(fh,dt,nh,ph,b0h)
        endif
c       do 12 i = 1, nh
c         write(*,*) ' i, ph = ', i, ph(i)
c12     continue

c---Get low pass filter poles if necessary

        if(nl.ne.0) then
          allocate(pl(nl)) 
          call lowpass(fl,dt,nl,pl,b0l)
        endif
c       do 14 i = 1, nh
c         write(*,*) ' i, pl = ', i, pl(i)
c14     continue

c-----Through with calculation of poles -
c  Copy over data for processing array
        do 15 i = 1, npts
          af(i) = of(i)
15      continue

c  now start filtering the data
c  Filter is implemented as a cascade of second order filters
c  high pass filter first use poles ph
c  Numerator polynomial is z**2 - 2*z + 1
c       write(*,*) ' nh, nl = ', nh, nl
        if (nh.ne.0) then
          do 10 i = 1, nh, 2
c  get first set of second order filter coefficients
c  from each pair of poles
            preal = dreal(ph(i))
            pimag = dimag(ph(i))
            a1 = -2.d0*preal
            a2 = preal*preal + pimag*pimag
            b1 = -2.d0
            b2 = 1.d0
            call filt (a1, a2, b1, b2, npts, af, af)
10        continue
c apply gain section
          do 20 i = 1, npts
            af(i) = b0h*af(i)
20        continue
        end if
c  apply low pass filter using poles pl
c  Numerator polynomial is z**2 + 2*z + 1

        if (nl.ne.0) then
          do 30 i = 1, nl, 2
            preal = dreal(pl(i))
            pimag = dimag(pl(i))
            a1 = -2.d0*preal
            a2 = preal*preal + pimag*pimag
            b1 = 2.d0
            b2 = 1.d0
            call filt (a1, a2, b1, b2, npts, af, af)
30        continue

c---apply gain section

          write(*,*) 'high gain is : ',b0h
          write(*,*) 'low gain is : ',b0l
          do 40 i = 1, npts
            af(i) = b0l*af(i)
40        continue
        end if
        if (allocated(pl)) deallocate(pl)
        if (allocated(ph)) deallocate(ph) 
        afin(1:npts) = sngl(af(1:npts)) !copy back 
        deallocate(af)
        deallocate(of) 
        return
        end

c**********************************************************************
c        filt (a1, a2, b1, b2, npts, fi, fo)
c**********************************************************************
c       Routine to apply a second order recursive filter to the data
c       denomonator polynomial is z**2 + a1*z + a2
c       numerator polynomial is z**2 + b1*z + b2
c           fi = input array
c           fo = output array
c           npts = number of points
c**********************************************************************


        subroutine filt (a1, a2, b1, b2, npts, fi, fo)
        implicit real*8 (a-h,o-z)
        real*8 fi(npts), fo(npts)
        integer npts
        real*8 d1, d2, out


c       double precision d1, d2, out
        integer i

        d1 = 0.d0
        d2 = 0.d0
        do 10 i = 1, npts
          out = fi(i) + d1
          d1 = b1*fi(i) - a1*out + d2
          d2 = b2*fi(i) - a2*out
          fo(i) = out
10      continue

        return
        end

c************************************************************************
c                   lowpass (fc,dt,n,p,b)
c************************************************************************
c
c    Routine to compute lowpass filter poles for Butterworth filter
c       fc = desired cutoff frequency
c       dt = sample rate in seconds
c       n = number of poles (MUST BE EVEN)
c       p = pole locations (RETURNED)
c       b = gain factor for filter (RETURNED)
c
c
c   Program calculates a continuous Butterworth low pass IIRs with required
c    cut off frequency.
c   This program is limited to using an even number of poles
c   Then a discrete filter is calculated utilizing the bilinear transform
c   Methods used here follow those in Digital Filters and Signal Processing
c   by Leland B. Jackson
c************************************************************************

        subroutine lowpass (fc,dt,n,p,b)
        implicit real*8 (a-h,o-z)
        complex*16 p(n)
        real*8 dt, b
        integer n

        integer i
        complex*16 one, x, y
        data pi /3.141592653589793d0/

c       Initialize variables
        wcp = 2.d0 * fc * pi
        wc = (2.d0/dt)*dtan(wcp*dt/2.d0)
        one = dcmplx (1.d0, 0.d0)

c       Calculate position of poles for continuous filter
        do 10 i = 1, n, 2
          preal = -wc*dcos(dfloat(i)*PI/dfloat(2*n))
          pimag =  wc*dsin(dfloat(i)*PI/dfloat(2*n))
          p(i) =   dcmplx(preal, pimag)
          p(i+1) = dcmplx(preal, -pimag)
10      continue

c       Calculate position of poles for discrete filter using
c       the bilinear transformation

        do 20 i = 1, n, 2
          p(i) = p(i)*dcmplx(dt/2.d0,0.d0)
          x = one + p(i)
          y = one - p(i)
          p(i) = x/y
          p(i+1) = dconjg(p(i))
20      continue

c       calculate filter gain

        b0 = 1.d0
        do 30 i = 1, n, 2
          x = one - p(i)
          y = one - p(i+1)
          x = x*y
          b0 = b0*4.d0/dreal(x)
30      continue
        b0 = 1.d0/b0
        b = b0

        return
        end


c************************************************************************
c                   highpass (fc,dt,n,p,b)
c************************************************************************
c
c    Routine to compute lowpass filter poles for Butterworth filter
c       fc = desired cutoff frequency
c       dt = sample rate in seconds
c       n = number of poles (MUST BE EVEN)
c       p = pole locations (RETURNED)
c       b = gain factor for filter (RETURNED)
c
c
c   Program calculates a continuous Butterworth highpass IIRs
c   First a low pass filter is calculated with required cut off frequency.
c   Then this filter is converted to a high pass filter
c   This program is limited to using an even number of poles
c   Then a discrete filter is calculated utilizing the bilinear transform
c   Methods used here follow those in Digital Filters and Signal Processing
c   by Leland B. Jackson
c
c************************************************************************

        subroutine highpass (fc,dt,n,p,b)
        implicit real*8 (a-h,o-z)
        complex*16    p(n)
        integer     n

        integer i
        complex*16 one, x, y

        data pi /3.141592653589793d0/
c       Initialize variables

        wcp = 2.d0 * fc * pi
        wc = (2.d0/dt)*dtan(wcp*dt/2.d0)
        alpha = dcos(wc*dt)
        one = dcmplx(1.d0,0.d0)

c       get poles for low pass filter

        call lowpass(fc,dt,n,p,b0)

c       now find poles for highpass filter

        do 10 i = 1, n, 2
          x = dcmplx(alpha,0.d0)*one - p(i)
          y = one - dcmplx(alpha,0.d0)*p(i)
          p(i) = x/y
          p(i+1) = dconjg(p(i))
10      continue

c      Calculate gain for high pass filter

        b0 = 1.d0
        do 20 i = 1, n, 2
          x = one + p(i)
          y = one + p(i+1)
          x = x*y
          b0 = b0*4.d0/dreal(x)
20      continue
        b0 = 1.d0/b0
        b = b0

        return
        end
