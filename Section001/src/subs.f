      

	subroutine zero(x,n1,n2)
c-----
c	zero the real array from elements n1 to n2
c-----
	real x(*)
	do 1000 i=n1,n2
		x(i) = 0.0
 1000	continue
	return
	end


       function ask(quest)
c
c   interactive i-o for real numbers
c
      character quest*(*)
      integer ounit
      character*8 myformat
      common /innout/ inunit,ounit
c
      ilen = len(quest)
      write(myformat,'(a2,i3.3,a3)')'(a',ilen,',$)'
c      
      write(ounit,myformat) quest
c     write(ounit,100) (quest(j:j),j=1,len(quest))
      read(inunit,*) anser
      ask=anser
      return
  100 format(80(a1,$))
      end
      subroutine asktxt(quest,answer)
c
c   interactive i-o for character strings
c      string returned may be a maximun of 64 characters
c
      character quest*(*)
      character*64 answer
      character*8 myformat
      integer ounit
      common /innout/ inunit,ounit
c
      ilen = len(quest)
      write(myformat,'(a2,i3.3,a3)')'(a',ilen,',$)'
c      
      write(ounit,myformat) quest
c
c     write(ounit,100) (quest(j:j),j=1,len(quest))
c
      read(inunit,200) answer
  100 format(80(a1,$))
  200 format(a64)
      return
      end
      integer function blank(file)
c
c   routine to find the position of the first blank (' ')
c   character in a character variable called file
c
c   len is a UNIX Ridge function to find the length of a character variable
c   then leng defines the maximum value of blank
c
      character*(*) file
      integer ounit
      common /innout/ inunit,ounit
      leng=len(file)
      do 1 i=1,leng
         if(file(i:i).ne.' ') go to 1
            blank=i-1
            return
    1 continue
c  
c   if no ' ' is found, set blank = leng   
c 
      blank=leng
      return
      end
      subroutine coef8(p,vp1,vs1,ro1,vp2,vs2,ro2,ncode,nd,rmod,rph)
c
c     the routine coef8 is designed for the computation of reflection
c     and transmission coefficients at a plane interface between two
c     homogeneous solid halfspaces or at a free surface of a homogeneous
c     solid halfspace.
c
c     the codes of individual coefficients are specified by the
c     following numbers
c     a/ interface between two solid halfspaces
c     p1p1...1       p1s1...2       p1p2...3       p1s2...4
c     s1p1...5       s1s1...6       s1p2...7       s1s2...8
c     b/ free surface (for ro2.lt.0.00001)
c     pp.....1       px.....5       px,pz...x- and z- components of the
c     ps.....2       pz.....6       coef.of conversion,incident p wave
c     sp.....3       sx.....7       sx,sz...x- and z- components of the
c     ss.....4       sz.....8       coef.of conversion,incident s wave
c
c     i n p u t   p a r a m e t e r s
c           p...ray parameter
c           vp1,vs1,ro1...parameters of the first halfspace
c           vp2,vs2,ro2...parameters of second halfspace. for the free
c                    surface take ro2.lt.0.00001,eg.ro2=0., and
c                    arbitrary values of vp2 and vs2
c           ncode...code of the computed coefficient
c           nd...=0  when the positive direction of the ray
c                    and the x-axis make an acute angle
c                =1  when the wave impinges on the interface
c                    against the positive direction of the x-axis
c
c     o u t p u t   p a r a m e t e r s
c           rmod,rph...modul and argument of the coefficient
c
c     n o t e s
c     1/ positive p...in the direction of propagation
c     2/ positive s...to the left from p
c     3/ time factor of incident wave ... exp(-i*omega*t)
c     4/ formulae are taken from cerveny ,molotkov, psencik, ray method
c        in seismology, pages 30-35. due to the note 2, the signs at
c        certain coefficients are opposite
c
c       written by v.cerveny,1976
c       modified by t.j. owens, 3/22/82 see comments in code
c
      complex b(4),rr,c1,c2,c3,c4,h1,h2,h3,h4,h5,h6,h,hb,hc
      dimension prmt(4),d(4),dd(4)
c
      if(ro2.lt.0.000001)go to 150
      prmt(1)=vp1
      prmt(2)=vs1
      prmt(3)=vp2
      prmt(4)=vs2
      a1=vp1*vs1
      a2=vp2*vs2
      a3=vp1*ro1
      a4=vp2*ro2
      a5=vs1*ro1
      a6=vs2*ro2
      q=2.*(a6*vs2-a5*vs1)
      pp=p*p
      qp=q*pp
      x=ro2-qp
      y=ro1+qp
      z=ro2-ro1-qp

      g1=a1*a2*pp*z*z
      g2=a2*x*x
      g3=a1*y*y
      g4=a4*a5
      g5=a3*a6
      g6=q*q*pp
      do 21 i=1,4
      dd(i)=p*prmt(i)
   21 d(i)=sqrt(abs(1.-dd(i)*dd(i)))

      if(dd(1).le.1..and.dd(2).le.1..and.dd(3).le.1..and.dd(4).le.1.)
     1go to 100
c
c     complex coefficients
      do 22 i=1,4
      if(dd(i).gt.1.)go to 23
      b(i)=cmplx(d(i),0.)
      go to 22
   23 b(i)= cmplx(0.,d(i))
   22 continue
      c1=b(1)*b(2)
      c2=b(3)*b(4)
      c3=b(1)*b(4)
      c4=b(2)*b(3)
      h1=g1
      h2=g2*c1
      h3=g3*c2
      h4=g4*c3
      h5=g5*c4
      h6=g6*c1*c2
      h=1./(h1+h2+h3+h4+h5+h6)
      hb=2.*h
      hc=hb*p
      go to (1,2,3,4,5,6,7,8),ncode
    1 rr=h*(h2+h4+h6-h1-h3-h5)
      go to 26
    2 rr=vp1*b(1)*hc*(q*y*c2+a2*x*z)
      if(nd.ne.0)rr=-rr
      go to 26
    3 rr=a3*b(1)*hb*(vs2*b(2)*x+vs1*b(4)*y)
      go to 26
    4 rr=-a3*b(1)*hc*(q*c4-vs1*vp2*z)
      if(nd.ne.0)rr=-rr
      go to 26
    5 rr=-vs1*b(2)*hc*(q*y*c2+a2*x*z)
      if(nd.ne.0)rr=-rr
      go to 26
    6 rr=h*(h2+h5+h6-h1-h3-h4)
      go to 26
    7 rr=a5*b(2)*hc*(q*c3-vp1*vs2*z)
      if(nd.ne.0)rr=-rr
      go to 26
    8 rr=a5*b(2)*hb*(vp1*b(3)*y+vp2*b(1)*x)
      go to 26
c     real coefficients
  100 e1=d(1)*d(2)
      e2=d(3)*d(4)
      e3=d(1)*d(4)
      e4=d(2)*d(3)
      s1=g1
      s2=g2*e1
      s3=g3*e2
      s4=g4*e3
      s5=g5*e4
      s6=g6*e1*e2
      s=1./(s1+s2+s3+s4+s5+s6)
      sb=2.*s
      sc=sb*p
      go to (101,102,103,104,105,106,107,108),ncode
  101 r=s*(s2+s4+s6-s1-s3-s5)
      go to 250
  102 r=vp1*d(1)*sc*(q*y*e2+a2*x*z)
      if(nd.ne.0)r=-r
      go to 250
  103 r=a3*d(1)*sb*(vs2*d(2)*x+vs1*d(4)*y)
      go to 250
  104 r=-a3*d(1)*sc*(q*e4-vs1*vp2*z)
      if(nd.ne.0)r=-r
      go to 250
  105 r=-vs1*d(2)*sc*(q*y*e2+a2*x*z)
      if(nd.ne.0)r=-r
      go to 250
  106 r=s*(s2+s5+s6-s1-s3-s4)
      go to 250
  107 r=a5*d(2)*sc*(q*e3-vp1*vs2*z)
      if(nd.ne.0)r=-r
      go to 250
  108 r=a5*d(2)*sb*(vp1*d(3)*y+vp2*d(1)*x)
      go to 250
c
c     earths surface,complex coefficients and coefficients of conversion
c
c   n o t e :
c
c    signs of coefficients at loops 162, 166, & 168 have been changed
c    from the originnal version of coef8 due to inconsistencies in
c    notation from the cerveny, et al book
c    3/22/82
c
  150 a1=vs1*p
      a2=a1*a1
      a3=2.*a2
      a4=2.*a1
      a5=a4+a4
      a6=1.-a3
      a7=2.*a6
      a8=2.*a3*vs1/vp1
      a9=a6*a6
      dd(1)=p*vp1
      dd(2)=p*vs1
      do 151 i=1,2
  151 d(i)=sqrt(abs(1.-dd(i)*dd(i)))
      if(dd(1).le.1..and.dd(2).le.1.)go to 200
      do 154 i=1,2
      if(dd(i).gt.1.)go to 155
      b(i)=cmplx(d(i),0.)
      go to 154
  155 b(i)= cmplx(0.,d(i))
  154 continue
      h1=b(1)*b(2)
      h2=h1*a8
      h=1./(a9+h2)
      go to (161,162,163,164,165,166,167,168),ncode
  161 rr=(-a9+h2)*h
      go to 26
  162 rr=-a5*b(1)*h*a6
      if(nd.ne.0)rr=-rr
      go to 26
  163 rr=a5*b(2)*h*a6*vs1/vp1
      if(nd.ne.0)rr=-rr
      go to 26
  164 rr=-(a9-h2)*h
      go to 26
  165 rr=a5*h1*h
      if(nd.ne.0)rr=-rr
      go to 26
  166 rr=-a7*b(1)*h
      go to 26
  167 rr=a7*b(2)*h
      go to 26
  168 rr=a5*vs1*h1*h/vp1
      if(nd.ne.0)rr=-rr
   26 z2=real(rr)
      z3=aimag(rr)
      if(z2.eq.0..and.z3.eq.0.)go to 157
      rmod=sqrt(z2*z2+z3*z3)
      rph=atan2(z3,z2)
      return
  157 rmod=0.
      rph=0.
      return
c
c     earths surface,real coefficients and coefficients of conversion
c   n o t e :
c
c    signs of coeficients at loops 202, 206, & 208 have been reversed
c    by t.j. owens because of inconsistencies w/sign conventions
c    3/22/82
c
  200 s1=d(1)*d(2)
      s2=a8*s1
      s=1./(a9+s2)
      go to (201,202,203,204,205,206,207,208),ncode
  201 r=(-a9+s2)*s
      go to 250
  202 r=-a5*d(1)*s*a6
      if(nd.ne.0)r=-r
      go to 250
  203 r=a5*d(2)*s*a6*vs1/vp1
      if(nd.ne.0)r=-r
      go to 250
  204 r=(s2-a9)*s
      go to 250
  205 r=a5*s1*s
      if(nd.ne.0)r=-r
      go to 250
  206 r=-a7*d(1)*s
      go to 250
  207 r=a7*d(2)*s
      go to 250
  208 r=a5*vs1*s1*s/vp1
      if(nd.ne.0)r=-r
  250 if(r.lt.0.)go to 251
      rmod=r
      rph=0.
      return
  251 rmod=-r
      rph=-3.14159
      return
      end
      subroutine coefsh(p,vs1,rho1,vs2,rho2,ncode,rmod,rph)
c
c  calculates sh-wave reflection and transmission coeficients
c    at a solid-solid interface of a free surface
c
c   solid-solid:  ncode =>
c
c    s1s1 = 1   s1s2 = 2
c
c   free surface: rho2 = 0.0 ncode =>
c
c    s1s1 = 1   free surface correction = 2
c
      complex p2,p4,d,h1,h2,rr
      if(rho2.lt..00001) go to 5
      a1=vs1*p
      a2=vs2*p
      b1=rho1*vs1
      b2=rho2*vs2
      g1=sqrt(abs(1. - a1*a1))
      g2=sqrt(abs(1. - a2*a2))
      p2=cmplx(g1,0.)
      p4=cmplx(g2,0.)
      if(a1.gt.1.) p2=cmplx(0.,g1)
      if(a2.gt.1.) p4=cmplx(0.,g2)
      h1=cmplx(b1,0.)*p2
      h2=cmplx(b2,0.)*p4
      d= h1 + h2
      go to (1,2) ncode
    1 rr=(h1-h2)/d
      go to 3
    2 rr=2.*h1/d
    3 z1=real(rr)
      z2=aimag(rr)
      if(z2.eq.0.) go to 4
      rmod=sqrt(z1*z1 + z2*z2)
      rph=atan2(z2,z1)
      return
    4 rmod=z1
      rph=0.
      return
c
c    free surface problem
c
    5 go to (6,7) ncode
    6 rmod=1.
      rph=0.
      return
    7 rmod=2.
      rph=0.
      return
      end
      subroutine coord(x,theta,delta,y,trans,same)
c
c  transforms a vector x in one coordinate system to a vector y
c    in another coordinate system defined by strike of theta
c    and a dip of delta where y' is the strike direction and
c    z' is the dip direction of the new system wrt to x and z of
c    the old system respectively
c
c    trans defines the direction of the transform --
c     if trans = 'local' then y will be in the primed system
c     if trans = 'globe' then y will be in the original system
c
c    if same = .true. then the transformation matrix is not recalculted
c                          from the previous call
c
      dimension x(3),y(3),a(3,3)
      character trans*5
      logical same
      integer ounit
      common /innout/ inunit,ounit
      common /cord/ a
      if(same) go to 4
      cost=cos(theta)
      sint=sin(theta)
      cosd=cos(delta)
      sind=sin(delta)
      a(1,1)=cost
      a(2,1)=-cosd*sint
      a(3,1)=sind*sint
      a(1,2)=sint
      a(2,2)=cosd*cost
      a(3,2)=-sind*cost
      a(1,3)=0.
      a(2,3)=sind
      a(3,3)=cosd
   4  if(trans.eq.'globe') go to 1
      if(trans.ne.'local') go to 5
      do 2 i=1,3
    2    y(i)=a(i,1)*x(1)+a(i,2)*x(2)+a(i,3)*x(3)
      return
    1 do 3 i=1,3
    3    y(i)=a(1,i)*x(1)+a(2,i)*x(2)+a(3,i)*x(3)
      return
    5 write(ounit,101) trans
  101 format(' trans = ',a5,' in coord, no transformation done')
      return
      end
      integer function daymo (dofy,month,day,year)
c
c function daymo determines the month and day
c of the month,given the year and day of year.
c it returns 1 if it was successful,0 otherwise.
c if dofy is not within legal limits,month and
c day will be returned as zero.
c
c
c calls:
c   lpyr
c
c      programmed by madeleine zirbes
c         september 15,1980
c
c day of year - input
      integer dofy
c month - output
      integer month
c day of month - output
      integer day
c year - input
      integer year
c
c day of year
      integer iday
c function
      integer lpyr
c number of days in month
      integer mdays(12)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31/
c
      iday = dofy
      if(.not.(iday.lt.1))goto 23000
         month = 0
         day = 0
         daymo = (0)
         return
c
23000 continue
      if(.not.(lpyr(year).eq.1))goto 23002
         mdays(2) = 29
         goto 23003
c     else
23002    continue
         mdays(2) = 28
23003 continue
c
      do 23004 month = 1,12
         day = iday
         iday = iday - mdays(month)
         if(.not.(iday.le.0))goto 23006
            daymo = (1)
            return
23006    continue
23004    continue
c
      month = 0
      day = 0
      daymo = (0)
      return
c
      end
      subroutine dfftr (x,nft,dirctn,delta)
c                                              a.shakal, 1/78, 15 jul 80
c           this subroutine does a fast fourier transform on a real
c        time series.  it requires 1/2 the storage and e1/2 the time
c        required by a complex fft.
c
c     forward transform, "call dfftr(x,nft,'forward',dt)":
c           input = x(1),x(2),..,x(nft) = real time series of nft points
c          output = x(1),x(2),..,x(nft+2) = nft/2+1 complex spectral poi
c        these spectral points are identical to the first nft/2+1 return
c        by subroutine fft (i.e., pos freq terms).  thus, the coefficien
c        at fj, the j-th frequency point (where fj = (j-1)*delf, j=1,nft
c        and delf = 1/(nft*dt)), is in x(i-1),x(i), where i=2j.  x(1) is
c        dc term, x(2) = 0 (because real time series), x(nft+1) is real
c        of nyquist coef, and x(nft+2) is imaginary part (0 because real
c        series).
c
c     inverse transform, "call dfftr(x,nft,'inverse',delf)":
c        input and output are interchanged.
c
c           if this subroutine is called with 'forward', and then with '
c        and delf of 1/(nft*dt), the original time series is recovered.
c        identical results (but for scaling) can be obtained by calling
c        fft(x,nft,isign), but in fft a real time series must be stored
c        complex array with zero imaginary parts, which requires 2*nft p
c        of array x.  also, the coefs returned by the fft will differ by
c        n-scaling, since fft's leave out the dt,delf of the approximate
c        integrations.  this subroutine calls fft.
c           this subroutine is a modification of the subroutine 'fftr',
c        written by c.frasier.  the principal modifications are:
c             1) the delt,delf of the integrations are included to make
c                a discrete approximation to the fourier transform.
c             2) the storage of the spectrum (on output if forward, or i
c                if inverse) has x(2) = zero, with the nyquist component
c                x(nft+1), with x(nft+2) = 0.
c
      logical forwrd, invrse
      character dirctn*7
      complex  csign, c1, c2, c3, speci, specj
      real x(nft+2)
      integer ounit
      common /innout/ inunit,ounit
      pi = 3.1415927
c
      call locast(dirctn,invrse,forwrd)
c
      nftby2 = nft/2
      if (.not.(forwrd)) go to 20001
c            forward transform..
      call fft (x,nftby2,-1)
      x1 = x(1)
      x(1) = x1 + x(2)
      x(2) = x1 - x(2)
      sign = -1.
      go to 20002
20001 if (.not.(invrse)) go to 10001
c            adjust nyquist element storage for inverse transform
      x(2) = x(nft+1)
      x(nft+1) = 0.
      sign = +1.
      go to 20002
10001 stop 'dirctn bad to dfftr'
c
c           manipulate elements as appropropriate for a 1/2 length
c        complex fft, after the forward fft, or before the inverse.
20002 piovrn = pi*sign/float(nftby2)
      csign = cmplx(0.,sign)
      do 10 i = 3,nftby2,2
      j = nft-i+2
      c1 = cmplx(x(i)+x(j), x(i+1)-x(j+1))
      c2 = cmplx(x(i)-x(j), x(i+1)+x(j+1))
      w = piovrn*float(i/2)
      c3 = cmplx(cos(w),sin(w))*c2
      speci = c1 + csign*c3
      x(i) = real(speci)/2.
      x(i+1) = aimag(speci)/2.
      specj = conjg(c1) + csign*conjg(c3)
      x(j) = real(specj)/2.
      x(j+1) = aimag(specj)/2.
   10 continue
      x(nftby2+2) = -x(nftby2+2)
      if (.not.(forwrd)) go to 20004
c            include dt of integration, for forward transform...
      dt = delta
      do 9000  i = 1,nft
 9000 x(i) = x(i)*dt
c            adjust storage of the nyquist component...
      x(nft+1) = x(2)
      x(nft+2) = 0.
      x(2) = 0.
      go to 20005
20004 if (.not.(invrse)) go to 10002
      x1 = x(1)
      x(1) = (x1+x(2))/2.
      x(2) = (x1-x(2))/2.
c            do the inverse transform...
      call fft (x,nftby2,+1)
c            in the inverse transform, include the df of the integration
c            and a factor of 2 because only doing half the integration
c            (i.e., just over the positive freqs).
      twodf = 2.*delta
      do 9002  i = 1,nft
 9002 x(i) = x(i)*twodf
10002 continue
20005 return
      end
      function diffr(x,y)
c
c   from lawson and hanson
c
      diffr=x-y
      return
      end
      function dot(x,y)
c
c  calculates the dot product of two vectors
c
      dimension x(1),y(1)
      z=0.
      do 1 i=1,3
    1    z=z + x(i)*y(i)
      dot=z
      return
      end
      integer function doy (month,day,year)
c
c function doy determines the day of the
c year,given the month,da ad year.
c if month or day are illegal,the return
c value of the function is zero.
c
c
c calls:
c   lpyr
c
c      programmed by madeleine zirbes
c         september 15,1980
c
c month - input
      integer month
c day of month - input
      integer day
c year - input
      integer year
c function
      integer lpyr
      integer inc
      integer ndays(12)
      data ndays /0,31,59,90,120,151,181,212,243,273,304,334/
c
      if(.not.(month .lt.1 .or. month .gt. 12))goto 23000
         doy = (0)
         return
23000 continue
      if(.not.(day .lt. 1 .or. day .gt. 31))goto 23002
         doy = (0)
         return
23002 continue
      if(.not.(lpyr(year).eq.1 .and. month .gt. 2))goto 23004
         inc = 1
         goto 23005
c     else
23004    continue
         inc = 0
23005 continue
      doy = (ndays(month) + day + inc)
      return
      end
      subroutine fft(data,nn,isign)
c                                              a.shakal, 1/78, 10 jul 80
c        cooley-tukey 'fast fourier trnasform' in ansi fortran 77.
c
c           transform(j) = sum {data(i)*w**u(i-1)*(j-1)e}, where i and
c        j run from 1 to nn, and w = exp(sign*twopi*sqrtu-1e/nn).
c        data is a one-dimensional complex array (i.e., the real and
c        imaginary parts of the data are located immediately adjacent
c        in storage, such as fortran places them) whose length nn is
c        a power of two.  isign is +1 or -1, giving the sign of the
c        transform.  transform values are returned in array data,
c        replacing the input data.  the time is proportional to
c        n*log2(n), rather than the non-fft n**2.  modified from the
c        fortran ii coding from n.brenner's mit-ll tech rept.
c
      real data(1)
      pi = 3.1415926
c
      n = 2*nn
      j = 1
      do 5 i = 1,n,2
      if (.not.(i .lt. j)) go to 20001
      tempr = data(j)
      tempi = data(j+1)
      data(j) = data(i)
      data(j+1) = data(i+1)
      data(i) = tempr
      data(i+1) = tempi
20001 m = n/2
    3 if (.not.(j .gt. m)) go to 20004
      j = j-m
      m = m/2
      if (m .ge. 2) go to 3
20004 j = j+m
   5  continue
c
c
      mmax = 2
    6 if (.not.(mmax .ge. n)) go to 20007
      return
20007 if (.not.(mmax .lt. n)) go to 10001
      istep = 2*mmax
      pibymx = pi*float(isign)/float(mmax)
c
      do 8 m = 1,mmax,2
      theta = pibymx*float(m-1)
      wr = cos(theta)
      wi = sin(theta)
      do 8 i = m,n,istep
      j = i + mmax
      tempr = wr*data(j) - wi*data(j+1)
      tempi = wr*data(j+1) + wi*data(j)
      data(j) = data(i) - tempr
      data(j+1) = data(i+1) - tempi
      data(i) = data(i) + tempr
      data(i+1) = data(i+1) + tempi
   8  continue
      mmax = istep
      go to 6
10001 continue
20008 return
      end
      function fsorce (isorfn,f,delay,kstrtd,a,b,t,wo,trap)
c
c           function to return the (complex) spectral value of a specifi
c        source function at frequency f, with a specifiable time delay.
c            # 1  s(t) = (a/dt)del(t)              s(w) = a
c            # 2  s(t) = triangle, a by t          s(w) = sqrd sinc(wt/4
c            # 3  s(t) = at*exp(-bt)               s(w) = a/(bb-ww,2wb)
c            # 4  s(t) = a*exp(-bt)*sin(ct)        s(w) = ac/(bb+cc-ww,2
c        calls to ask, sinc.                 a.shakal 10/76, 10/5/78, 8
c
      complex fsorce
      dimension trap(4)
      integer ounit
      common /innout/ inunit,ounit
      data e,twopi/2.7182818,6.2831853/
c
c
      w = twopi*f
      tdelay = delay
      if(kstrtd .gt. 0) go to 20
c
c        initialization of constants.
c
      go to (11,12,13,14,15,15,16), isorfn
      stop ' isorfn ng in fsorce'
c
c           source function # 1.  s(w)=flat, s(t)=spike.  amp(t)=amp(w)/
c
 11   a = ask ('s(w) = a, s(t) = spike of amp a/dt. a= ')
      go to 20
c
c           source function # 2.  triangle of width t, amplitude a.  spe
c        is s(w) = (at/2)*sinc(wt/4)**2*exp(-iwt/2).  it peaks at w = 0,
c        where s(w) = at/2.
c
   12 write(ounit,110)
  110 format(1x,'s(t)=triangle of width t. s(w)=sinc sqrd')
      t= ask('width, t (sec)= ')
      atmp = ask ('amp (in time or freq, +/-) = ')
      if (atmp .gt. 0.) a = atmp*t/2.
      if (atmp .lt. 0.) a = -atmp
      go to 20
c
c           source function # 3.  non-oscillating pulse, bell-shaped spe
c        centered at 0.  s(t)=atexp(-bt), s(w)=a/(b+iw)**2.  max amp in
c        is amp in freq*b/e.  max amp in time is at t =  1/b, at which s
c        a/(b*e).  max in freq is at w = 0, where s(w) = a/b**2.
c
   13 write(ounit,111)
  111 format(1x,'non-oscillating pulse. bell-shaped spectrum',
     *          ' centered at 0')
      trise = ask ('rise time= ')
      b = 1./trise
      atmp = ask ('amp (in time or freq, +/-) = ')
      if (atmp .gt. 0.) a = b*e*atmp
      if (atmp .lt. 0.) a = -b*b*atmp
      go to 20
c
c           source function # 4.  damped sinusoid.  s(t) = aexp(-bt)sin(
c
 14   fo = ask ('damped sinusoid. freq fo = ')
      wo = twopi*fo
      b = ask ('decay, if other than wo/5 ')
      if (b .eq. 0.) b = wo/5.
      atmp = ask ('amp (in time or freq, +/-) = ')
      if (atmp .gt. 0.) tmax = atan(wo/b)/wo
      if (atmp .gt. 0.) a = atmp/(exp(-b*tmax)*sin(wo*tmax))
      if (atmp .lt. 0.) a = -2.*b*atmp
      go to 20
   15 trap(1)=ask('trapezoid.  rise time= ')
      trap(2)=ask('falloff time= ')
      trap(3)=ask('total duration= ')
      trap(4)=ask('amplitude= ')
      go to 20
c
c   source function #4  gaussian
c
   16 a=ask('gaussian function exp(-w*w/4*a*a) a = ')
c
c
 20   go to (21,22,23,24,25,25,26), isorfn
 21   fsorce = cmplx(a, 0.)
      go to 30
 22   if(w.gt.1.e-06) go to 225
      fsorce=cmplx(a,0.)
      tdelay=tdelay+t/2.
      go to 30
 225  fsorce = cmplx(a*(sinc(w*t/4.))**2,0.)
      tdelay = tdelay + t/2.
      go to 30
 23   fsorce = cmplx(a, 0.)/cmplx(b**2-w**2, 2.*w*b)
      go to 30
 24   fsorce = cmplx(a*wo, 0.)/cmplx(b**2+wo**2-w**2, 2.*b*w)
      go to 30
   25 if(w.gt.1e-06) go to 31
      fsorce=cmplx((trap(4)/2.)*(2.*trap(3)-trap(2)-trap(1)),0.)
      go to 30
   26 gauss=-w*w/(4.*a*a)
      fsorce=cmplx(exp(gauss),0.)
      go to 30
   31 ral=(cos(w*trap(1))-1.)/trap(1) +
     *    (cos(w*(trap(3)-trap(2)))-cos(w*trap(3)))
     *    /trap(2)
      aim=sin(w*trap(1))/trap(1) +
     *    (sin(w*(trap(3)-trap(2)))-sin(w*trap(3)))
     *    /trap(2)
      fac=trap(4)/(w*w)
      fsorce=cmplx(fac*ral,-fac*aim)
c
c        apply a time delay by shift theorem
c
 30   if(tdelay .ne. 0.) fsorce = fsorce*exp(cmplx(0., -w*tdelay))
      kstrtd = 1
      return
      end
      subroutine g1(a,b,cos,sin,sig)
c
c   from lawson and hanson
c
      zero=0.
      one=1.
      if(abs(a).le.abs(b)) go to 10
      xr=b/a
      yr=sqrt(one+xr**2)
      cos=sign(one/yr,a)
      sin=cos*xr
      sig=abs(a)*yr
      return
   10 if(b) 20,30,20
   20 xr=a/b
      yr=sqrt(one+xr**2)
      sin=sign(one/yr,b)
      cos=sin*xr
      sig=abs(b)*yr
      return
   30 sig=zero
      cos= zero
      sin=one
      return
      end
      subroutine g2(cos,sin,x,y)
c
c   from lawson and hanson
c
      xr=cos*x+sin*y
      y=-sin*x+cos*y
      x=xr
      return
      end
      subroutine h12(mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)
c
c     from lawson and hanson
c
      dimension u(iue,m),c(1)
      double precision sm,b
      one=1.
      if(0.ge.lpivot.or.lpivot.ge.l1.or.l1.gt.m) return
      cl=abs(u(1,lpivot))
      if(mode.eq.2) go to 60
         do 10 j=l1,m
   10    cl=amax1(abs(u(1,j)),cl)
      if(cl) 130,130,20
   20 clinv=one/cl
      sm=(dble(u(1,lpivot))*clinv)**2
         do 30 j=l1,m
   30    sm=sm+(dble(u(1,j))*clinv)**2
      sm1=sm
      cl=cl*sqrt(sm1)
      if(u(1,lpivot)) 50,50,40
   40 cl=-cl
   50 up=u(1,lpivot)-cl
      u(1,lpivot)=cl
      go to 70
   60 if(cl) 130,130,70
   70 if(ncv.le.0) return
      b=dble(up)*u(1,lpivot)
      if(b) 80,130,130
   80 b=one/b
      i2=1-icv+ice*(lpivot-1)
      incr=ice*(l1-lpivot)
         do 120 j=1,ncv
         i2=i2+icv
         i3=i2+incr
         i4=i3
         sm=c(i2)*dble(up)
            do 90 i=l1,m
            sm=sm+c(i3)*dble(u(1,i))
   90       i3=i3+ice
         if(sm) 100,120,100
  100    sm=sm*b
         c(i2)=c(i2)+sm*dble(up)
            do 110 i=l1,m
            c(i4)=c(i4)+sm*dble(u(1,i))
  110       i4=i4+ice
  120    continue
  130 return
      end
      function iask(quest)
c
c      interactive i-o for integers
c
      integer answer,ounit
      character quest*(*)
      common /innout/ inunit,ounit
      write(ounit,100) (quest(j:j),j=1,len(quest))
      read(inunit,*) answer
      iask=answer
      return
  100 format(80(a1,$))
      end
      subroutine iniocm
c
c     subroutine to initialize all variables in /tjocm/
c
c **********************************************************************
c
c common block info for link with subroutine sacio
c
      real instr
      integer year,jday,hour,min,isec,msec,idf
      character*8 sta,cmpnm,evnm,kdf
      common /tjocm/ dmin,dmax,dmean,year,jday,hour,min,isec,msec,sta,
     *            cmpnm,az,cinc,evnm,baz,delta,rayp,depth,decon,agauss,
     *              c,tq,instr,dlen,begin,t0,t1,t2
c
c **************************************************************************
c
      data df,idf,kdf,zero/-12345.,-12345,'        ',0./
      dmin=zero
      dmax=zero
      dmean=zero
      year=idf
      jday=idf
      hour=idf
      min=idf
      isec=idf
      msec=idf
      sta=kdf
      cmpnm=kdf
      az=df
      cinc=df
      evnm=kdf
      baz=df
      delta=df
      rayp=df
      depth=df
      decon=zero
      agauss=df
      c=df
      tq=df
      instr=df
      dlen=df
      begin=df
      t0=df
      t1=df
      t2=df
      return
      end
      subroutine juli(yr,jday,month,iday,monum)
c
c     converts julian day to month dat,year
c       month in a3 format
c       iday = day of month
c       monum = number of month
c
      dimension mon(12),num(12)
      integer yr
      character mon*3,month*3
      data mon/'jan','feb','mar','apr','may','jun','jul','aug',
     *         'sep','oct','nov','dec'/,
     *     num/31,28,31,30,31,30,31,31,30,31,30,31/
      iday=jday
      ind=0
      do 1 i=1,12
         ind=ind + num(i)
         if(yr.ne.(yr/4)*4) go to 2
            if(i.ne.2) go to 2
              ind=ind+1
    2    if(iday.gt.ind) go to 1
            month=mon(i)
            iday=num(i)-(ind-iday)
            monum = i
            go to 3
    1 continue
    3 return
      end
      subroutine locast(dirctn,invrse,forwrd)
      character dirctn*7
      logical forwrd,invrse
      integer ounit
      common /innout/ inunit,ounit
      if(dirctn.eq.'forward') go to 1
      if(dirctn.eq.'inverse') go to 2
      write(ounit,100)dirctn
  100 format(1x,a7,2x,'is meaningless to dfftr, use forward or inverse
     *only')
      invrse=.false.
      forwrd=.false.
      return
    1 invrse=.false.
      forwrd=.true.
      return
    2 invrse=.true.
      forwrd=.false.
      return
      end
      integer function lpyr(year)
c
c function lpyr determines if year
c is a leap year.
c
c this function uses the intrinsic
c function mod. if your machine
c does not supply this function,
c make one -
c mod(i,j) = iabs(i - (i/j)*j)
c
c
c calls:
c   mod - intrinsic funtion
c
c      programmed by madeleine zirbes
c         september 15,1980
c
c year - input
      integer year
      if(.not.(mod(year,400).eq.0))goto 23000
         lpyr = (1)
         return
23000 continue
      if(.not.(mod(year,4) .ne. 0))goto 23002
         lpyr = (0)
         return
23002 continue
      if(.not.(mod(year,100).eq.0))goto 23004
         lpyr = (0)
         return
23004 continue
      lpyr = (1)
      return
      end
      subroutine max(x,n,xmax)
      dimension x(1)
      xmax=0.
      do 1 i=1,n
         if(abs(x(i)).gt.xmax) xmax=abs(x(i))
    1 continue
      return
      end
      subroutine minmax(x,npts,min,max,mean)
      dimension x(1)
      real min,max,mean
      min=9.0e+19
      max=-9.0e+19
      mean=0.
      do 1 i=1,npts
           if(x(i).gt.max) max=x(i)
           if(x(i).lt.min) min=x(i)
           mean=mean + x(i)
    1 continue
      mean=mean/float(npts)
      return
      end
      function npowr2(n)
c
c finds the next power of 2 .ge.n
c
      ipowr=alog10(2.*float(n)-1.)/.301029996
      if(n.eq.1) ipowr=1
      npowr2=2**ipowr
      return
      end
      subroutine qrbd(ipass,q,e,nn,v,mdv,nrv,c,mdc,ncc)
c
c    from lwason and hanson
c
      logical wntv,havers,fail
      dimension q(nn),e(nn),v(mdv,nn),c(mdc,ncc)
      zero=0.
      one=1.
      two=2.
      n=nn
      ipass=1
      if(n.le.0) return
      n10=n*10
      wntv=nrv.gt.0
      havers=ncc.gt.0
      fail=.false.
      nqrs=0
      e(1)=zero
      dnorm=zero
         do 10 j=1,n
   10    dnorm=amax1(abs(q(j))+abs(e(j)),dnorm)
         do 200 kk=1,n
         k=n+1-kk
   20    if(k.eq.1) go to 50
         if(diffr(dnorm+q(k),dnorm)) 50,25,50
   25    cs=zero
         sn=-one
            do 40 ii=2,k
            i=k+1-ii
            f=-sn*e(i+1)
            e(i+1)=cs*e(i+1)
            call g1(q(i),f,cs,sn,q(i))
            if(.not.wntv) go to 40
               do 30 j=1,nrv
   30          call g2(cs,sn,v(j,i),v(j,k))
   40       continue
   50       do 60 ll=1,k
            l=k+1-ll
            if(diffr(dnorm+e(l),dnorm)) 55,100,55
   55       if(diffr(dnorm+q(l-1),dnorm)) 60,70,60
   60       continue
         go to 100
   70    cs=zero
         sn=-one
            do 90 i=l,k
            f=-sn*e(i)
            e(i)=cs*e(i)
            if(diffr(dnorm+f,dnorm)) 75,100,75
   75       call g1(q(i),f,cs,sn,q(i))
            if(.not.havers) go to 90
               do 80 j=1,ncc
   80          call g2(cs,sn,c(i,j),c(l-1,j))
   90       continue
  100    z=q(k)
         if(l.eq.k) go to 170
         x=q(l)
         y=q(k-1)
         g=e(k-1)
         h=e(k)
         f=((y-z)*(y+z)+(g-h)*(g+h))/(two*h*y)
         g=sqrt(one+f**2)
         if(f.lt.zero) go to 110
         t=f+g
         go to 120
  110    t=f-g
  120    f=((x-z)*(x+z)+h*(y/t-h))/x
         cs=one
         sn=one
         lp1=l+1
            do 160 i=lp1,k
            g=e(i)
            y=q(i)
            h=sn*g
            g=cs*g
            call g1(f,h,cs,sn,e(i-1))
            f=x*cs+g*sn
            g=-x*sn+g*cs
            h=y*sn
            y=y*cs
            if(.not.wntv) go to 140
               do 130 j=1,nrv
  130          call g2(cs,sn,v(j,i-1),v(j,i))
  140       call g1(f,h,cs,sn,q(i-1))
            f=cs*g+sn*y
            x=-sn*g+cs*y
            if(.not.havers) go to 160
               do 150 j=1,ncc
  150          call g2(cs,sn,c(i-1,j),c(i,j))
  160       continue
         e(l)=zero
         e(k)=f
         q(k)=x
         nqrs=nqrs+1
         if(nqrs.le.n10) go to 20
         fail=.true.
  170    if(z.ge.zero) go to 190
         q(k)=-z
         if(.not.wntv) go to 190
            do 180 j=1,nrv
  180       v(j,k)=-v(j,k)
  190    continue
  200    continue
      if(n.eq.1) return
         do 210 i=2,n
         if(q(i).gt.q(i-1)) go to 220
  210    continue
      if(fail) ipass=2
      return
  220    do 270 i=2,n
         t=q(i-1)
         k=i-1
            do 230 j=i,n
            if(t.ge.q(j)) go to 230
            t=q(j)
            k=j
  230       continue
         if(k.eq.i-1) go to 270
         q(k)=q(i-1)
         q(i-1)=t
         if(.not.havers) go to 250
            do 240 j=1,ncc
            t=c(i-1,j)
            c(i-1,j)=c(k,j)
  240       c(k,j)=t
  250    if(.not.wntv) go to 270
            do 260 j=1,nrv
            t=v(j,i-1)
            v(j,i-1)=v(j,k)
  260       v(j,k)=t
  270    continue
      if(fail) ipass=2
      return
      end
      subroutine rdlyrs (sitefl,nlyrs,title,vp,vs,rho,h,
     *qp,qs,strike,dip,iflag,ier)
c
c        rdlyr - read a layered-medium file
c     arguments:
c        iflag = -1 a dipping layer model is read in with
c                   the strike and dip of the bottom of the ith layer
c                   given in the input file
c              = 0  an infinite q, flat layered model is assumed
c              = 1  a finite q, flat-layered model is assumed, qp & qs
c                   must be non-zero in the file
c        ier = 0 unless the layered numbers are screwed up
c
      character*32 sitefl,title
      logical yes,yesno
      real  vp(1),vs(1),qp(1),qs(1),rho(1),h(1),strike(1),dip(1)
      integer ounit
      common /innout/ inunit,ounit
      iu=20
      yes=yesno('List the site model? (y or n) ')
      open(unit=iu,file=sitefl)
      rewind iu
      read (iu,100) nlyrs,title
      ier=0
      write(ounit,104) sitefl,title,nlyrs
      if(yes) write(ounit,105)
      do 6 i=1,nlyrs
      read(iu,101) k,vp(i),vs(i),rho(i),h(i),qpk,qsk,theta,delta
      if(k.eq.i) go to 4
        write(ounit,102)
        ier=1
        close(unit=iu)
        return
    4 if(iflag) 1,2,3
    1 strike(i)=theta
      dip(i)=delta
      if(i.ne.nlyrs) go to 5
         if((theta.ne.0)) write(ounit,103)
         if((delta.ne.0)) write(ounit,103)
         go to 5
    2 qp(i)=-1.
      qs(i)=-1.
      go to 5
    3 qp(i)=qpk
      qs(i)=qsk
      if(qp(i).lt.0) qp(i)=.75*(vp(i)/vs(i))**2*qs(i)
    5 if(yes) write(ounit,106) i,vp(i),vs(i),rho(i),h(i),
     *                         qpk,qsk,theta,delta
    6 continue
      close(unit=iu)
      return
  100 format(i3,a10)
  101 format(i3,1x,8f8.2)
  102 format(' layers out of order in rdlyrs *******')
  103 format(' warning -- strike and dip of half space were given')
  104 format(' file: ',a10,' model: ',a10,2x,i3,' layers ')
  105 format(/,/,' lyr     vp      vs     rho      h     qp',
     *           '      qs     strike    dip')
  106 format(1x,i3,1x,8f8.2)
      end
      subroutine rotate(x,m,n,baz,az,npts)
c
c   rotates horz. components of x into radial & tangential components
c      given the back azimuth and their orientations
c
c   baz   = back azimuth from station to source in degrees
c   az(i) = + direction of each horz. comp.
c
c      conventions --
c
c          radial is positive away form source
c          tangential is positive clockwise from + radial direction
c
      dimension x(m,n),az(3)
      integer ounit
      common /innout/ inunit,ounit
      rad(deg)=deg/57.295779
      azck=az(2) + 90.
      diff=abs(azck - az(3))
      if(diff.lt..01) go to 1
      if(diff.lt.179.) write(ounit,100) az(2),az(3),azck,diff
         do 2 i=1,npts
    2     x(i,3)=-x(i,3)
    1 a=sin(rad(baz) - rad(az(2)))
      b=cos(rad(baz) - rad(az(2)))
      do 4 i=1,npts
         radial = -x(i,2)*b - x(i,3)*a
         trans  =  x(i,2)*a - x(i,3)*b
         x(i,2) = radial
         x(i,3) = trans
    4 continue
      return
  100 format(' p r o b l e m   i n   r o t a t e ',/,
     *       1x,'+ horz direction: 1= ',f8.4,' 2= ',f8.4,/,
     *       1x,'for proper rotation 2= ',f8.4,' difference = ',f8.4)
      end
      subroutine sacio(file,x,np,dt,inout)
      include './mach'
      include './hdr'
      integer year,jday,hour,min,isec,msec,ounit
      character*8 sta,compnm,evnm
      character file*64
      dimension x(1)
      common /tjocm/ dmin,dmax,dmean,year,jday,hour,min,isec,msec,sta,
     *              compnm,caz,cinc,evnm,bz,del,rayp,dep,decon,agauss,
     *              c,tq,rinstr,dlen,btimey,ty0,ty1,ty2
      common /innout/ inunit,ounit
c **********************************************************************
c
c common block info for link with subroutine sacio
c
c
c **************************************************************************
c
c   parameters are:
c
c      dmin,dmax,dmean = min,max, and mean of data read or written
c      year,jday,hour,min,isec,msec = gmt reference time (all integers)
c      sta = station name (a8)
c      compnm = component name (a8)
c      caz = orientation of component (wrt north)
c      cinc = inclination of component (wrt vertical)
c      evnm = name of event (a8)
c      bz  = back azimuth of from station to event
c      del = distance from station to event in degrees
c      rayp  = ray parameter of arriving phase in j-b earth
c      dep = depth of event in kilometers
c      decon = if = 1. indicates data has been source equalized
c      agauss = width of gaussian used in source equalization if decon = 1.
c      c     = trough filler used in source equalization if decon = 1.
c      tq = t/q value used in synthetic (if used)
c      rinstr = 1. if response of 15-100 system has been put into synthetic
c      dlen  = length of data in secs.
c      btimey = time 0f 1st data point wrt gmt reference time (in secs)
c      ty0,ty1,ty2 = user defined times wrt gmt reference (in secs)
c
c ****************************************************************************
c
c    call to sacio is:
c                      call sacio(file,x,np,dt,inout)
c
c    where file = file to be read or written
c          x    = data array to be used
c          np   = number of points in x
c          dt   = sampling rate for x
c          inout = +1 for reading a sac file
c                = -1 for writing a sac file
c
c ******************************************************************************
      if(inout.lt.0) goto 1
c
c  read a sac file
c
c     call inicm
      call inihdr()
      call rsac1(file,x,nnp,bb,ddt,100000,nerr)
      if(nerr.eq.0) go to 2
      write(ounit,100) nerr
  100 format(' nerr = ',i2,' in sacio read')
      return
C    2 sta=kstnm
    2	call getkhv('KSTNM'  , sta, ierr)
	call getnhv('NPTS'   , np,  ierr)
	b = bb
	dt = ddt
	np =nnp
C      np=npts
C      dt=delta
	call getfhv('DELTA', dt, ierr)
C      dmin=depmin
C      dmax=depmax
C      dmean=depmen
	call getfhv('DEPMIN', dmin, ierr)
	call getfhv('DEPMAX', dmax, ierr)
	call getfhv('DEPMEN', demean, ierr)
C      year=nzyear
C      hour=nzhour
C      jday=nzjday
C      min=nzmin
C      isec=nzsec
C      msec=nzmsec
	call getnhv('NZYEAR', year, ierr)
	call getnhv('NZHOUR', hour, ierr)
	call getnhv('NZJDAY', jday, ierr)
	call getnhv('NZMIN', min, ierr)
	call getnhv('NZSEC', isec, ierr)
	call getnhv('NZMSEC', msec, ierr)
	call getfhv('B',b,ierr)
	call getfhv('E',e,ierr)
      btimey=b
      dlen=e-b
	call getkhv('KCMPNM', compnm, ierr)
C      compnm=kcmpnm
C      caz=cmpaz
C      cinc=cmpinc
	call getfhv('CMPAZ   ', caz, ierr)
	call getfhv('CMPINC  ', cinv, ierr)
	call getkhv('KEVNM', evnm, ierr)
C      evnm=kevnm
C      bz=baz
C      del=gcarc
	call getfhv('GCARC', del, ierr)
	call getfhv('BAZ', bz, ierr)
C      rayp=user0
C      dep=user1
C      agauss=user2
C      c=user3
	call getfhv('USER0   ', rayp, ierr)
	call getfhv('USER1   ', dep, ierr)
	call getfhv('USER2   ', agauss, ierr)
	call getfhv('USER3   ', c, ierr)
      if(c.gt.0.) decon=1.0
C      tq=user4
C      rinstr=user5
	call getfhv('USER4   ', tq, ierr)
	call getfhv('USER5   ', rinstr, ierr)
C      ty0=t0
C      ty1=t1
C      ty2=t2
	call getfhv('T0      ', ty0, ierr)
	call getfhv('T1      ', ty1, ierr)
	call getfhv('T2      ', ty2, ierr)
      return
    1 call inihdr()
c
c write a sac file
c
      call newhdr()
      if(year.lt.1960) go to 4
      if(decon.lt..0001) go to 5
      b = 0.0
      b=btimey
      go to 6
    5 nzyear=year
      nzhour=hour
      nzjday=jday
      nzmin=min
      nzsec=isec
      nzmsec=msec
    4 b=btimey
    6 if(abs(ty0).gt..00001) t0=ty0
      if(abs(ty1).gt..00001) t1=ty1
      if(abs(ty2).gt..00001) t2=ty2
      kcmpnm=compnm
      cmpaz=caz
      cmpinc=cinc
      kevnm=evnm
      baz=bz
      if(del.gt..0001) gcarc=del
      user0=rayp
      if(dep.gt..0001) user1=dep
      if(decon.lt..0001) go to 3
      user2=agauss
      user3=c
    3 if(tq.gt..0001) user4=tq
      if(rinstr.gt..0001) user5=rinstr
    	call setkhv('KSTNM'  , sta, ierr)
	call setnhv('NPTS'   , np,  ierr)
	call setfhv('DELTA', dt, ierr)
	call setfhv('DEPMIN', dmin, ierr)
	call setfhv('DEPMAX', dmax, ierr)
	call setfhv('DEPMEN', dmean, ierr)
	call setnhv('NZYEAR', year, ierr)
	call setnhv('NZHOUR', nzhour, ierr)
	call setnhv('NZJDAY', nzjday, ierr)
	call setnhv('NZMIN', nzmin, ierr)
	call setnhv('NZSEC', nzsec, ierr)
	call setnhv('NZMSEC',nzmsec, ierr)
	call setfhv('B',b,ierr)
	call setfhv('E',e,ierr)
	call setkhv('KCMPNM', kcmpnm, ierr)
	call setfhv('CMPAZ   ', cmpaz, ierr)
	call setfhv('CMPINC  ', cmpinc, ierr)
	call setkhv('KEVNM', evnm, ierr)
	call setfhv('GCARC', del, ierr)
	call setfhv('BAZ', baz, ierr)
	call setfhv('USER0   ', user0, ierr)
	call setfhv('USER1   ', user1, ierr)
	call setfhv('USER2   ', user2, ierr)
	call setfhv('USER3   ', user3, ierr)
	call setfhv('USER4   ', user4, ierr)
	call setfhv('USER5   ', user5, ierr)
	call setfhv('T0      ', t0, ierr)
	call setfhv('T1      ', t1, ierr)
	call setfhv('T2      ', t2, ierr)
	call setlhv('LEVEN   ',.true.,ierr)
      call wsac0(file,xxx,x,nerr)
C      call wsac1(file,x,np,b,dt,nerr)
      if(nerr.eq.0) return
      write(ounit,101) nerr
  101 format(' nerr = ',i2,' in sacio write')
      return
      end
      subroutine seisio(freq,peak,xr,xi,inout)
c     15-100 system
c     wwssn instrument constants from u. chandra, bssa 1970 vol 60
c     pp 539-563
c     peak magnifications are 350,700,1400,2800,5600
c
c     if inout = +1 response is put in
c     if inout = -1 response is removed
c
      if(freq.lt.0.005) freq = 0.005
    8 we = 6.2831853*freq
      index = (peak + 1.)/375.
      go to (1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5),index
    1 fmag = 278.
      sigma = 0.003
      go to 6
    2 fmag = 556.
      sigma = 0.013
      go to 6
    3 fmag = 1110.
      sigma = 0.047
      go to 6
    4 fmag = 2190.
      sigma = 0.204
      go to 6
    5 fmag = 3950.
      sigma = 0.805
    6 zeta = 0.93
      zeta1=1.
      wn=.418879
      wn1 = .062831853
      ar= (we*we-wn*wn)*(we*we-wn1*wn1)-4.*zeta*zeta1*wn*wn1*(1.-sigma)
     1*we*we
             ai=2.*we*(zeta1*wn1*(wn*wn-we*we)+zeta*wn*(wn1*wn1-we*we))
      if(inout.eq.+1) go to 7
      factor = 1./(fmag*we*we*we)
      xr = - ai * factor
      xi = ar * factor
      return
    7 factor=fmag*we*we*we/(ai*ai + ar*ar)
      xr=-factor*ai
      xi=-factor*ar
      return
      end
      function sinc(x)
      integer ounit
      common /innout/ inunit,ounit
      if(abs(x).le.1.0e-08) go to 1
      sinc=sin(x)/x
      return
    1 sinc=0.
      write(ounit,100)
  100 format(' arg sinc = 0, sinc set to 0')
      return
      end
      subroutine sit2(psvsh,c,freq,resp,alfa,beta,qp,qs,rho,thik,nlyrs)
c
c           compute the site response for nlyrs-1 dissipative layers ove
c        a halfspace for a p, sv or sh wave incident from the halfspace.
c        frequency domain solution, yeilding the complex response at a
c        given frequency and phase veloctiy.
c        this solution sub-divides layers for which p or q are too
c        large at this freq into equal sub-layers, to prevent machine ov
c        flow problems.
c
c          arguments...
c        psvsh = 1,2,3 for an incident p, sv or sh wave.
c
c        freq,c - prescribed freq (hz) & horizontal phase velocity (c is
c            not restricted to be greater than alfa or beta)
c
c        resp - the complex response in the u,v,w directions.  the respo
c            is the ratio of the free surface to incident displacements
c            velocities, etc), the 'crustal transfer functions' of haske
c            (1962).  for p or sv solution, resp(1) & (3) contain the ho
c            zontal and vertical (u and w) response.  for sh, resp(2)
c            contains the  horizontal (v) response.
c
c        alfa,beta,qp,qs,rho and thik contain the medium properties for
c            layers 1 thru nlyrs (the halfspace).  thik(nlyrs) not refrn
c            input qs(i) as negative if want no dissipation in that laye
c            if qp(i) is input negative, it is calculated from qs(i),
c            assuming qkappa is infinite (eg., stacey,1969).
c
c        note: this is the same solution method as subroutine 'crust', b
c        made specific for the crustal transfer function solution.
c         original routine by:  a. shakal  9/78
c         modified by t.j. owens 7/81
c
      logical psvwav,shwave,test
      integer psvsh,hsize,hafsiz
      real k,alfa(1),beta(1),qp(1),qs(1),rho(1),thik(1)
c
c        complex declarations are in alphabetical order...
c
      complex cosp,cospq1,cosq,d,f,fg,g,gama,gamas,gamas2,g
     &amas3,   gm1,gsm1,gm1gs1,h(4,4),hn(4,2),hprod(4,2)
     &,             i,j(4,2),          kalfa,kbeta,nu,nualfa,nubeta,nupr
     &od,one,p,q,   resp(3),ro,s,shterm,sinp,sinq,two,zero
      integer ounit
      common /innout/ inunit,ounit
      data twopi,eps/6.2831853,.001/,  expmax/30./,  i,zero/(0.,1.),(0.,
     &0.)/,   one,two/(1.,0.),(2.,0.)/
      nu(c,vel) = cmplx(sqrt(abs((c/vel)**2 -1.)), 0.)
c
c
      w = twopi*freq
      if(freq .eq. 0.) w = 1.0e-6
      signw = sign(1.,w)
      k = w/c
      shwave = psvsh .eq. 3
      psvwav = psvsh .le. 2
      hsize = 4
      if(shwave) hsize = 2
      hafsiz = hsize/2
      do 9000  ii=1,3
 9000 resp(ii) = zero
      p = zero
c
c        check that none of the q's are input as zero  (cc pull for spee
c
      do 10 n = 1,nlyrs
      if(abs(qp(n)) .gt. eps  .and.  abs(qs(n)) .gt. eps) go to 10
      write (ounit,108) n
      stop
  10  continue
  108     format('site: q=0 in layr',i3)
c
c
c           calculate hn, the matrix of propagation down to the top of
c        the halfspace:  hn = h(nlyrs-1)*h(nlyrs-2)*...*h(2)*h(1)
c
c
      do 40 n = 1,nlyrs
c
c        values depending on s-wave velocity...
c
      qsinv = amax1(0., 1./qs(n))
      if(qs(n) .lt. 0.) kbeta = cmplx(k, 0.)
      if(qs(n) .gt. 0.) kbeta = cmplx(k, 0.)/csqrt(cmplx(1., signw*qsinv
     &))
      nubeta = nu(c,beta(n))
      if(c .lt. beta(n)) nubeta = -i*nubeta
      gama = cmplx(2.*(beta(n)/c)**2, 0.)
      ro = cmplx(rho(n), 0.)
      if(shwave) shterm = two*kbeta/(ro*gama*nubeta)
      if(shwave) go to 22
c
c        values depending on p-wave velocity...
c
      qpinv = 1./qp(n)
      if(qp(n) .lt. 0.) qpinv = 1.33*(beta(n)/alfa(n))**2*qsinv
      if(qpinv .le. 0.) kalfa = cmplx(k, 0.)
      if(qpinv .gt. 0.) kalfa = cmplx(k, 0.)/csqrt(cmplx(1., signw*qpinv
     &))
      nualfa = nu(c,alfa(n))
      if(c .lt. alfa(n)) nualfa = -i*nualfa
      nuprod = nualfa*nubeta
      s = kalfa/kbeta
      gamas = gama*s
      gamas2 = gamas*gamas
      gamas3 = gamas2*gamas
      gm1 = gama - one
      gsm1 = gamas*s - one
      gm1gs1 = gm1*gsm1
      f = one + gamas - gamas*s
      g = one - gama + gamas
      fg = f*g
 22   if(n .eq. nlyrs .and. nlyrs .gt. 1) go to 42
c
c        if p or q will cause sinh or cosh to be too large for the machi
c        given this layer thickness, divide layer into 'nparts' equal pa
c
      q = kbeta*nubeta*cmplx(thik(n), 0.)
      if(psvwav) p = kalfa*nualfa*cmplx(thik(n), 0.)
      nparts = amax1(abs(aimag(p)),abs(aimag(q)))/expmax +1.
      if(nparts.gt.1) write(ounit,985) w,nparts,n
  985 format(1x,e15.6,2i8)
      if(nparts .gt. 1) q = q/cmplx(float(nparts), 0.)
      if(nparts .gt. 1 .and. psvwav) p = p/cmplx(float(nparts), 0.)
      sinq =csin(q)
      cosq =ccos(q)
      if(shwave) go to 26
      sinp = csin(p)
      cosp = ccos(p)
      cospq1 = cosp - cosq
c
c        compute h, the 4x4 transfer matrix of this layer, analogous
c        to a(m) of haskell(1953).    for p-sv problem...
c
      h(1,1) = (gamas*cosp - gsm1*cosq)/f
      h(2,1) = -i*(gamas*nualfa*sinp + gsm1*sinq/nubeta)/f
      h(3,1) = ro*gama*gsm1*(cospq1)/(kbeta*f)
      h(4,1) = i*ro*(nualfa*gamas2*sinp + gm1gs1*sinq/nubeta)/(kbeta
     &*f)
      h(1,2) = i*(gm1*sinp/nualfa + gamas*nubeta*sinq)/g
      h(2,2) = -(gm1*cosp - gamas*cosq)/g
      h(3,2) = i*ro*(gm1gs1*sinp/nualfa + gamas2*nubeta*sinq)/(kalfa
     &*g)
      h(4,2) = ro*gamas*gm1*(cospq1)/(kbeta*g)
      if(n .eq. 1 .and. nparts .eq. 1) go to 30
      h(1,3) = -kalfa*(cospq1)/(ro*f)
      h(2,3) = i*kalfa*(nualfa*sinp + sinq/nubeta)/(ro*f)
      h(3,3) = -(gsm1*cosp - gamas*cosq)/f
      h(4,3) = -i*s*(gamas*nualfa*sinp + gm1*sinq/nubeta)/f
      h(1,4) = i*kbeta*(sinp/nualfa + nubeta*sinq)/(ro*g)
      h(2,4) = -kbeta*(cospq1)/(ro*g)
      h(3,4) = i*(gsm1*sinp/nualfa + gamas*nubeta*sinq)/(s*g)
      h(4,4) = (gamas*cosp - gm1*cosq)/g
      go to 30
c
c        for sh problem...
c
 26   h(1,1) = cosq
      h(2,1) = i*sinq/shterm
      h(1,2) = i*shterm*sinq
      h(2,2) = h(1,1)
c
c        multiply to obtain the matrix of propagation down thru layer n.
c        only need 1st col of hn for sh, and the 1st 2 cols, with the 1s
c        col of hhn, for psv.  if nparts > 1, do multiplication for each
c        of the nparts this layer is divided into.
c
   30 if(n .gt. 1) go to 34
      do 9004  jj=1,hafsiz
      do 9004  ii=1,hsize
 9004 hn(ii,jj) = h(ii,jj)
      if(nlyrs .eq. 1) go to 42
      nparts = nparts - 1
      if(nparts .eq. 0) go to 40
 34   do 36 npart = 1,nparts
      do 9010  jj = 1,hafsiz
      do 9010  ii = 1,hsize
 9010 hprod(ii,jj) = zero
      do 9012  jj = 1,hafsiz
      do 9012    ii = 1,hsize
      do 9012  kk = 1,hsize
 9012 hprod(ii,jj) = hprod(ii,jj) + h(ii,kk)*hn(kk,jj)
      do 9014  jj=1,hafsiz
      do 9014  ii=1,hsize
 9014 hn(ii,jj) = hprod(ii,jj)
 36   continue
 40   continue
c
c
c           solve for the 'crustal transfer functions', the ratio of the
c        free surface displacements (or velocities, etc) to the incident
c        (i.e., fullspace) displacements (or veloc.).  using left half
c        of j = einvrs*hn.
c
c
 42   if(shwave) go to 50
      do 48 m = 1,2
      j(3,m) = (gsm1*hn(1,m)/kbeta - s*hn(3,m)/ro)/(nubeta*f)
      j(4,m) = (gamas*s*hn(2,m)/kalfa + hn(4,m)/ro)/g
      j(1,m) = (gama*hn(1,m)/kbeta - hn(3,m)/ro)/f
 48   j(2,m) = -(gm1*hn(2,m)/kalfa + hn(4,m)/(ro*s))/(nualfa*g)
      d = (j(1,1)-j(2,1))*(j(3,2)-j(4,2))
     *   +(j(1,2)-j(2,2))*(j(4,1)-j(3,1))
      if(psvsh.eq.1) go to 51
c
c      incident sv wave.  calc usv, wsv of haskell's notation.
c
      resp(1) = two*(j(1,2) -j(2,2))/(kbeta*nubeta*d)
      resp(3) = two*(j(2,1) -j(1,1))/(kbeta*d)
      return
c
c      incident p wave.  calculate up, wp of haskell's notation.
c
   51 resp(1)=two*(j(3,2)-j(4,2))/(kalfa*d)
      resp(3)=two*(j(3,1)-j(4,1))/(kalfa*nualfa*d)
      return
c
c        incident sh wave.  calculate vsh.
c
 50   resp(2) = two/(hn(1,1) +shterm*hn(2,1))
      return
      end
      subroutine svdrs(a,mda,mm,nn,b,mdb,nb,s)
c
c
c    lawson and Hanson singular value decomposition routine
c
c	s occupies 3*n cells
c	a occupies m*n cells
c	b occupies m*nb cells
c
c	
      dimension a(mda,nn),b(mdb,nb),s(nn,3)
      integer ounit
      common /innout/ inunit,ounit
      zero=0.
      one=1.
      n=nn
      if(n.le.0.or.mm.le.0) return
      j=n
   10 continue
         do 20 i=1,mm
         if(a(i,j)) 50,20,50
   20    continue
      if(j.eq.n) go to 40
         do 30 i=1,mm
   30    a(i,j)=a(i,n)
   40 continue
      a(1,n)=j
      n=n-1
   50 continue
      j=j-1
      if(j.ge.1) go to 10
      ns=0
      if(n.eq.0) go to 240
      i=1
      m=mm
   60 if(i.gt.n.or.i.ge.m) go to 150
      if(a(i,i)) 90,70,90
   70    do 80 j=1,n
         if(a(i,j)) 90,80,90
   80    continue
      go to 100
   90 i=i+1
      go to 60
  100 if(nb.le.0) go to 115
         do 110 j=1,nb
         t=b(i,j)
         b(i,j)=b(m,j)
  110    b(m,j)=t
  115    do 120 j=1,n
  120    a(i,j)=a(m,j)
      if(m.gt.n) go to 140
         do 130 j=1,n
  130    a(m,j)=zero
  140 continue
      m=m-1
      go to 60
  150 continue
c
c   end .. special for zero rows and columns
c   begin .. svd alogritm
c
      l=min0(m,n)
         do 170 j=1,l
         if(j.ge.m) go to 160
         call h12(1,j,j+1,m,a(1,j),1,t,a(1,j+1),1,mda,n-j)
         call h12(2,j,j+1,m,a(1,j),1,t,b,1,mdb,nb)
  160    if(j.ge.n-1) go to 170
         call h12(1,j+1,j+2,n,a(j,1),mda,s(j,3),a(j+1,1),mda,1,m-j)
  170    continue
      if(n.eq.1) go to 190
         do 180 j=2,n
         s(j,1)=a(j,j)
  180    s(j,2)=a(j-1,j)
  190 s(1,1)=a(1,1)
      ns=n
      if(m.ge.n) go to 200
      ns=m+1
      s(ns,1)=zero
      s(ns,2)=a(m,m+1)
  200 continue
         do 230 k=1,n
         i=n+1-k
         if(i.ge.n-1) go to 210
         call h12(2,i+1,i+2,n,a(i,1),mda,s(i,3),a(1,i+1),1,mda,n-i)
  210       do 220 j=1,n
  220       a(i,j)=zero
  230    a(i,i)=one
      call qrbd(ipass,s(1,1),s(1,2),ns,a,mda,n,b,mdb,nb)
      go to (240,310), ipass
  240 continue
      if(ns.ge.n) go to 260
      nsp1=ns+1
         do 250 j=nsp1,n
  250     s(j,1)=zero
  260 continue
      if(n.eq.nn) return
      np1=n+1
         do 280 j=np1,nn
         s(j,1)=a(1,j)
            do 270 i=1,n
  270       a(i,j)=zero
  280    continue
         do 300 k=np1,nn
         i=s(k,1)
         s(k,1)=zero
            do 290 j=1,nn
            a(k,j)=a(i,j)
  290       a(i,j)=zero
         a(i,k)=one
  300    continue
      return
  310 write(ounit,320)
      stop
  320 format(49h convergence failure in qr bidiagonal svd routine)
      end
      logical function yesno(quest)
c
c   interactive i-o for logical variables
c    yesno must be declared logical in calling program
c
      character quest*(*),answer*1
      logical lanswr
      integer ounit
      character*8 myformat
      common /innout/ inunit,ounit
c
      ilen = len(quest)
      write(myformat,'(a2,i3.3,a3)')'(a',ilen,',$)'
c      
      write(ounit,myformat) quest
c     write(ounit,100) (quest(j:j),j=1,len(quest))

      read(inunit,200) answer
      lanswr=.false.
      if(answer.eq.'y') lanswr=.true.
      yesno=lanswr
c 100 format(80(a1,$))
  100 format(72(a1))
  200 format(a1)
      return
      end
