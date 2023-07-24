c-----------------------------------------------------
c     FFT subroutine
c-------------------------------------------------------   
      subroutine fourt(data,nn,ndim,isign,iform,work)
c-----------------------------------------------------------------------
c     data     --- the array used to hold the real and imaginary parts
c                  of the data on input 
c                  and the transform values on output.  
c     nn(ndim) --- integer array containing the dimensions of data
c     isign    --- -1 for a forward transform (exponential sign is -) 
c                  +1 for an inverse transform (sign is +).  
c     iform    --- 1 if the data and the transform values are complex;
c                  0 if the data are real but the transform values 
c                    are complex. If it is 0, the imaginary parts of 
c                    the data should be set to zero.  
c     work     --- an array used for working storage.  
c                  0 if all the dimensions of the data are powers of 
c		   two. 
c-----------------------------------------------------------------------
c
c     the cooley-tukey fast
c     transform(j1,j2,...) = sum(data(i1,i2,...)*w1**((i1-1)*(j1-1))
c                                 *w2**((i2-1)*(j2-1))*...),
c     sqrt(-1)/nn(1)), etc.  there is no limit on the dimensionality
c     (number of subscripts) of the data array.  if an inverse
c     where i1 and j1 run from 1 to nn(1) and w1=exp(isign*2*pi*
c     transform (isign=+1) is performed upon an array of transformed
c     (isign=-1) data, the original data will reappear,
c     multiplied by nn(1)*nn(2)*...  the array of input data may be
c     real or complex, at the programmers option, with a saving of
c     up to forty per cent in running time for real over complex.
c     (for fastest transform of real data, nn(1) should be even.)
c     the transform values are always complex, and are returned in the
c     original array of data, replacing the input data.  the length
c     of each dimension of the data array may be any integer.  the
c     program runs faster on composite integers than on primes, and is
c     particularly fast on numbers rich in factors of two.
c     timing is in fact given by the following formula.  let ntot be the
c     total number of points (real or complex) in the data array, that
c     is, ntot=nn(1)*nn(2)*...  decompose ntot into its prime factors,
c     such as 2**k2 * 3**k3 * 5**k5 * ...  let sum2 be the sum of all
c     the factors of two in ntot, that is, sum2 = 2*k2.  let sumf be
c     the sum of all other factors of ntot, that is, sumf = 3*k3+5*k5+..
c     the time taken by a multidimensional transform on these ntot data
c     is t = t0 + ntot*(t1+t2*sum2+t3*sumf).  on the cdc 3300 (floating
c     point add time = six microseconds), t = 3000 + ntot*(600+40*sum2+
c     175*sumf) microseconds on complex data.
c
c     implementation of the definition by summation will run in a time
c     proportional to ntot*(nn(1)+nn(2)+...).  for highly composite ntot
c     the savings offered by this program can be dramatic.  a one-dimen-
c     sional array 4000 in length will be transformed in 4000*(600+
c     40*(2+2+2+2+2)+175*(5+5+5)) = 14.5 seconds versus about 4000*
c     4000*175 = 2800 seconds for the straightforward technique.
c
c     the calling sequence is--
c     call fourt(data,nn,ndim,isign,iform,work)
c
c     data is the array used to hold the real and imaginary parts
c     of the data on input and the transform values on output.  it
c     is a multidimensional floating point array, with the real and
c     imaginary parts of a datum stored immediately adjacent in storage
c     (such as fortran iv places them).  the extent of each dimension
c     is given in the integer array nn, of length ndim.  isign is -1
c     to indicate a forward transform (exponential sign is -) and +1
c     for an inverse transform (sign is +).  iform is +1 if the data and
c     the transform values are complex.  it is 0 if the data are real
c     but the transform values are complex.  if it is 0, the imaginary
c     parts of the data should be set to zero.  as explained above, the
c     transform values are always complex and are stored in array data.
c     work is an array used for working storage.  it is not necessary
c     if all the dimensions of the data are powers of two.  in this case
c     it may be replaced by 0 in the calling sequence.  thus, use of
c     powers of two can free a good deal of storage.  if any dimension
c     is not a power of two, this array must be supplied.  it is
c     floating point, one dimensional of length equal to twice the
c     largest array dimension (i.e., nn(i) ) that is not a power of
c     two.  therefore, in one dimension for a non power of two,
c     work occupies as many storage locations as data.  if supplied,
c     work must not be the same array as data.  all subscripts of all
c     arrays begin at one.
c
c     the fast
c     nature of the data beyond the usual restriction that
c     the data form one cycle of a periodic function.  they are--
c     1.  the number of input data and the number of transform values
c     must be the same.
c     2. considering the data to be in the time domain,
c     they must be equi-spaced at intervals of dt.  further, the trans-
c     form values, considered to be in frequency space, will be equi-
c     spaced from 0 to 2*pi*(nn(i)-1)/(nn(i)*dt) at intervals of
c     2*pi/(nn(i)*dt) for each dimension of length nn(i).  of course,
c     dt need not be the same for every dimension.
c
c     there are no error messages or error halts in this program. the
c     program returns immediately if ndim or any nn(i) is less than one.
c     example 1.  three-dimensional forward
c     complex array dimensioned 32 by 25 by 13 in fortran iv.
c     dimension data(32,25,13),work(50),nn(3)
c     complex data
c     data nn/32,25,13/
c     do 1 i=1,32
c     do 1 j=1,25
c     do 1 k=1,13
c  1  data(i,j,k)=complex value
c     call
c
c     example 2.  one-dimensional forward transform of a real array of
c     length 64 in fortran ii.
c     dimension data(2,64)
c     do 2 i=1,64
c     data(1,i)=real part
c  2  data(2,i)=0.
c     call
c
c     program by norman brenner from the basic program by charles
c     rader (both of mit lincoln laboratory).  may 1967.  the idea
c     for the digit reversal was suggested by ralph alter (also mit ll).
c     this is the fastest and most versatile version of the fft known
c     to the author.  a program called
c     performs the fast
c     fortran.  it is about one third as long and restricts the
c     dimensions of the input array (which must be complex) to be powers
c     of two.  another program, called
c     runs two thirds as fast on a one-dimensional complex array whose
c     length is a power of two.
c
c     reference-- 
c     fast
c     g. sande, 1966 fall joint computer conference.
c
c     the work reported in this document was performed at lincoln lab-
c     oratory, a center for research operated by massachusetts institute
c     of technology, with the support of the u.s. air force under
c     contract af 19(628)-5167.
c
c     the fast
c
      dimension data(*), nn(ndim),ifact(32),work(*)
	pi=4.*atan(1.)
	twopi=pi+pi
c      twopi=6.283185307
      rthlf=.70710 67812
      np0=0
      nprev=0
C      if(ndim-1)920,1,1
      if((ndim-1).lt.0) go to 920
1     ntot=2
      do 2 idim=1,ndim
C      if(nn(idim))920,920,2
      if(nn(idim).le.0) go to 920
2     ntot=ntot*nn(idim)
c
c     main loop for each dimension
c
      np1=2
      do 910 idim=1,ndim
      n=nn(idim)
      np2=np1*n
C      if(n-1)920,900,5
      if((n-1).le.0) go to 920
c
c     is n a power of two and if not, what are its factors
c
5     m=n
      ntwo=np1
      if=1
      idiv=2
10    iquot=m/idiv
      irem=m-idiv*iquot
C      if(iquot-idiv)50,11,11
      if((iquot-idiv).lt.0) go to 50
C11    if(irem)20,12,20
11    if(irem.ne.0) go to 20
12    ntwo=ntwo+ntwo
      ifact(if)=idiv
      if=if+1
      m=iquot
      go to 10
20    idiv=3
      inon2=if
30    iquot=m/idiv
      irem=m-idiv*iquot
C      if(iquot-idiv)60,31,31
      if((iquot-idiv).lt.0) go to 60
C31    if(irem)40,32,40
31    if(irem.ne.0) go to 40
32    ifact(if)=idiv
      if=if+1
      m=iquot
      go to 30
40    idiv=idiv+2
      go to 30
50    inon2=if
C      if(irem)60,51,60
      if(irem.ne.0) go to 60
51    ntwo=ntwo+ntwo
      go to 70
60    ifact(if)=m
70    non2p=np2/ntwo
c
c     separate
c        1. complex transform
c        2. real transform for the 2nd, 3rd, etc. dimension.  method--
c           transform half the data, supplying the other half by con-
c           jugate symmetry.
c        3. real transform for the 1st dimension, n odd.  method--
c           set the imaginary parts to zero.
c        4. real transform for the 1st dimension, n even.  method--
c           transform a complex array of length n/2 whose real parts
c           are the even numbered real values and whose imaginary parts
c           are the odd numbered real values.  separate and supply
c           the second half by conjugate symmetry.
c
      ifmin = 1
      i1rng = np1
      if (iform .le. 0 .and. idim .lt. 4) go to 71
      icase = 1
      go to 100
c
71    if (idim .le. 1) go to 72
      icase = 2
      i1rng = np0 * (1 + nprev / 2)
      go to 100
c
72    if (ntwo .gt. np1) go to 73
      icase = 3
      go to 100
c
73    icase = 4
      ifmin = 2
      ntwo = ntwo / 2
      n = n / 2
      np2 = np2 / 2
      ntot = ntot / 2
c
      i = - 1
      do 80 j = 1, ntot
         i = i + 2
         data(j) = data(i)
80    continue
c
c     shuffle data by bit reversal, since n=2**k.  as the shuffling
c     can be done by simple interchange, no working array is needed
c
C100   if(non2p-1)101,101,200
100   if((non2p-1).gt.0) go to 200
101   np2hf=np2/2
      j=1
      do 150 i2=1,np2,np1
C      if(j-i2)121,130,130
      if((j-i2).ge.0) go to 130
121   i1max=i2+np1-2
      do 125 i1=i2,i1max,2
      do 125 i3=i1,ntot,np2
      j3=j+i3-i2
      tempr=data(i3)
      tempi=data(i3+1)
      data(i3)=data(j3)
      data(i3+1)=data(j3+1)
      data(j3)=tempr
125   data(j3+1)=tempi
130   m=np2hf
C140   if(j-m)150,150,141
140   if((j-m).le.0) go to 150
141   j=j-m
      m=m/2
C      if(m-np1)150,140,140
      if((m-np1).ge.0) go to 140
150   j=j+m
      go to 300
c
c     shuffle data by digit reversal for general n
c
200   nwork=2*n
      do 270 i1=1,np1,2
      do 270 i3=i1,ntot,np2
      j=i3
      do 260 i=1,nwork,2
C      if(icase-3)210,220,210
      if((icase-3).ne.0) go to 220
210   work(i)=data(j)
      work(i+1)=data(j+1)
      go to 240
220   work(i)=data(j)
      work(i+1)=0.
240   ifp2=np2
      if=ifmin
250   ifp1=ifp2/ifact(if)
      j=j+ifp1
C      if(j-i3-ifp2)260,255,255
      if((j-i3-ifp2).lt.0) go to 260
255   j=j-ifp2
      ifp2=ifp1
      if=if+1
C      if(ifp2-np1)260,260,250
      if((ifp2-np1).ge.0) go to 250
260   continue
      i2max=i3+np2-np1
      i=1
      do 270 i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
270   i=i+2
c
c     main loop for factors of two.
c     w=exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)
c     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2).
c
C300   if(ntwo-np1)600,600,305
300   if((ntwo-np1).le.0) go to 600
305   np1tw=np1+np1
      ipar=ntwo/np1
C310   if(ipar-2)350,330,320
310   if((ipar-2).lt.0) then
        go to 350
      elseif(ipar.eq.2) then
        go to 330
      else
      endif
320   ipar=ipar/4
      go to 310
330   do 340 i1=1,i1rng,2
      do 340 k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=data(k2)
      tempi=data(k2+1)
      data(k2)=data(k1)-tempr
      data(k2+1)=data(k1+1)-tempi
      data(k1)=data(k1)+tempr
340   data(k1+1)=data(k1+1)+tempi
350   mmax=np1
C360   if(mmax-ntwo/2)370,600,600
360   if((mmax-ntwo/2).ge.0) go to 600
370   lmax=max0(np1tw,mmax/2)
      do 570 l=np1,lmax,np1tw
      m=l
C      if(mmax-np1)420,420,380
      if((mmax-np1).le.0) go to 420
380   theta=-twopi*float(l)/float(4*mmax)
C      if(isign)400,390,390
      if(isign.lt.0) go to 400
390   theta=-theta
400   wr=cos(theta)
      wi=sin(theta)
410   w2r=wr*wr-wi*wi
      w2i=2.*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,i1rng,2
      kmin=i1+ipar*m
C      if(mmax-np1)430,430,440
      if((mmax-np1).gt.0) go to 440
430   kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
C      if(kstep-ntwo)460,460,530
      if((kstep-ntwo).gt.0) go to 530
460   do 520 k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
C      if(mmax-np1)470,470,480
      if((mmax-np1).gt.0) go to 480
470   u1r=data(k1)+data(k2)
      u1i=data(k1+1)+data(k2+1)
      u2r=data(k3)+data(k4)
      u2i=data(k3+1)+data(k4+1)
      u3r=data(k1)-data(k2)
      u3i=data(k1+1)-data(k2+1)
C      if(isign)471,472,472
      if(isign.ge.0) go to 472
471   u4r=data(k3+1)-data(k4+1)
      u4i=data(k4)-data(k3)
      go to 510
472   u4r=data(k4+1)-data(k3+1)
      u4i=data(k3)-data(k4)
      go to 510
480   t2r=w2r*data(k2)-w2i*data(k2+1)
      t2i=w2r*data(k2+1)+w2i*data(k2)
      t3r=wr*data(k3)-wi*data(k3+1)
      t3i=wr*data(k3+1)+wi*data(k3)
      t4r=w3r*data(k4)-w3i*data(k4+1)
      t4i=w3r*data(k4+1)+w3i*data(k4)
      u1r=data(k1)+t2r
      u1i=data(k1+1)+t2i
      u2r=t3r+t4r
      u2i=t3i+t4i
      u3r=data(k1)-t2r
      u3i=data(k1+1)-t2i
C      if(isign)490,500,500
      if(isign.ge.0) go to 500
490   u4r=t3i-t4i
      u4i=t4r-t3r
      go to 510
500   u4r=t4i-t3i
      u4i=t3r-t4r
510   data(k1)=u1r+u2r
      data(k1+1)=u1i+u2i
      data(k2)=u3r+u4r
      data(k2+1)=u3i+u4i
      data(k3)=u1r-u2r
      data(k3+1)=u1i-u2i
      data(k4)=u3r-u4r
520   data(k4+1)=u3i-u4i
      kdif=kstep
      kmin=4*(kmin-i1)+i1
      go to 450
530   continue
      m=m+lmax
C      if(m-mmax)540,540,570
      if((m-mmax).gt.0) go to 570
C540   if(isign)550,560,560
540   if(isign.ge.0) go to 560
550   tempr=wr
      wr=(wr+wi)*rthlf
      wi=(wi-tempr)*rthlf
      go to 410
560   tempr=wr
      wr=(wr-wi)*rthlf
      wi=(tempr+wi)*rthlf
      go to 410
570   continue
      ipar=3-ipar
      mmax=mmax+mmax
      go to 360
c
c     main loop for factors not equal to two.
c     w=exp(isign*2*pi*sqrt(-1)*(j1+j2-i3-1)/ifp2)
c
C600   if(non2p-1)700,700,601
600   if((non2p-1).le.0) go to 700
601   ifp1=ntwo
      if=inon2
610   ifp2=ifact(if)*ifp1
      theta=-twopi/float(ifact(if))
C      if(isign)612,611,611
      if(isign.lt.0) go to 612
611   theta=-theta
612   wstpr=cos(theta)
      wstpi=sin(theta)
      do 650 j1=1,ifp1,np1
      thetm=-twopi*float(j1-1)/float(ifp2)
C      if(isign)614,613,613
      if(isign.lt.0) go to 614
613   thetm=-thetm
614   wminr=cos(thetm)
      wmini=sin(thetm)
      i1max=j1+i1rng-2
      do 650 i1=j1,i1max,2
      do 650 i3=i1,ntot,np2
      i=1
      wr=wminr
      wi=wmini
      j2max=i3+ifp2-ifp1
      do 640 j2=i3,j2max,ifp1
      twowr=wr+wr
      j3max=j2+np2-ifp2
      do 630 j3=j2,j3max,ifp2
      jmin=j3-j2+i3
      j=jmin+ifp2-ifp1
      sr=data(j)
      si=data(j+1)
      oldsr=0.
      oldsi=0.
      j=j-ifp1
620   stmpr=sr
      stmpi=si
      sr=twowr*sr-oldsr+data(j)
      si=twowr*si-oldsi+data(j+1)
      oldsr=stmpr
      oldsi=stmpi
      j=j-ifp1
C      if(j-jmin)621,621,620
      if((j-jmin).gt.0) go to 620
621   work(i)=wr*sr-wi*si-oldsr+data(j)
      work(i+1)=wi*sr+wr*si-oldsi+data(j+1)
630   i=i+2
      wtemp=wr*wstpi
      wr=wr*wstpr-wi*wstpi
640   wi=wi*wstpr+wtemp
      i=1
      do 650 j2=i3,j2max,ifp1
      j3max=j2+np2-ifp2
      do 650 j3=j2,j3max,ifp2
      data(j3)=work(i)
      data(j3+1)=work(i+1)
650   i=i+2
      if=if+1
      ifp1=ifp2
C      if(ifp1-np2)610,700,700
      if((ifp1-np2).lt.0) go to 610
c
c     complete a real transform in the 1st dimension, n even, by con-
c     jugate symmetries.
c
700   go to (900,800,900,701),icase
701   nhalf=n
      n=n+n
      theta=-twopi/float(n)
C      if(isign)703,702,702
      if(isign.lt.0) go to 703
702   theta=-theta
703   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      imin=3
      jmin=2*nhalf-1
      go to 725
710   j=jmin
      do 720 i=imin,ntot,np2
      sumr=(data(i)+data(j))/2.
      sumi=(data(i+1)+data(j+1))/2.
      difr=(data(i)-data(j))/2.
      difi=(data(i+1)-data(j+1))/2.
      tempr=wr*sumi+wi*difr
      tempi=wi*sumi-wr*difr
      data(i)=sumr+tempr
      data(i+1)=difi+tempi
      data(j)=sumr-tempr
      data(j+1)=-difi+tempi
720   j=j+np2
      imin=imin+2
      jmin=jmin-2
      wtemp=wr*wstpi
      wr=wr*wstpr-wi*wstpi
      wi=wi*wstpr+wtemp
C725   if(imin-jmin)710,730,740
725   if((imin-jmin).lt.0) then
         go to 710
      elseif((imin-jmin).gt.0) then
         go to 740
      else
      endif
C730   if(isign)731,740,740
730   if(isign.ge.0) go to 740
731   do 735 i=imin,ntot,np2
735   data(i+1)=-data(i+1)
740   np2=np2+np2
      ntot=ntot+ntot
      j=ntot+1
      imax=ntot/2+1
745   imin=imax-2*nhalf
      i=imin
      go to 755
750   data(j)=data(i)
      data(j+1)=-data(i+1)
755   i=i+2
      j=j-2
C	if(i-imax)750,760,760
	if((i-imax).lt.0) go to 750
760   data(j)=data(imin)-data(imin+1)
      data(j+1)=0.
C      if(i-j)770,780,780
      if((i-j).lt.0) then
        go to 770
      else
        go to 780
      endif
765   data(j)=data(i)
      data(j+1)=data(i+1)
770   i=i-2
      j=j-2
C      if(i-imin)775,775,765
      if((i-imin).gt.0) go to 765
775   data(j)=data(imin)+data(imin+1)
      data(j+1)=0.
      imax=imin
      go to 745
780   data(1)=data(1)+data(2)
      data(2)=0.
      go to 900
c
c     complete a real transform for the 2nd, 3rd, etc. dimension by
c     conjugate symmetries.
c
C800   if(i1rng-np1)805,900,900
800   if((i1rng-np1).ge.0) go to 900
805   do 860 i3=1,ntot,np2
      i2max=i3+np2-np1
      do 860 i2=i3,i2max,np1
      imax=i2+np1-2
      imin=i2+i1rng
      jmax=2*i3+np1-imin
C      if(i2-i3)820,820,810
      if((i2-i3).le.0) go to 820
810   jmax=jmax+np2
C820   if(idim-2)850,850,830
820   if((idim-2).le.0) go to 850
830   j=jmax+np0
      do 840 i=imin,imax,2
      data(i)=data(j)
      data(i+1)=-data(j+1)
840   j=j-2
850   j=jmax
      do 860 i=imin,imax,np0
      data(i)=data(j)
      data(i+1)=-data(j+1)
860   j=j-np0
c
c     end of loop on each dimension
c
900   np0=np1
      np1=np2
910   nprev=n
920   return
      end


