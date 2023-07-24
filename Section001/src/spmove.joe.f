      program psmout
c
c  *************************************************************************
c
      parameter(MAXPOINTS=16388, MAXPOINTS2=MAXPOINTS*2)
      dimension d2(MAXPOINTS),caz(3),ampud(23),ampew(23),ampsn(23)
      character eqfile*64,outfil*64,comp(3,2)*6,outc(2)*5,knm(2)*8
      complex data(MAXPOINTS,3),zero1
      double precision gnorm
      logical yes,yesno,rel
      integer blank,ounit
c **********************************************************************
c
c common block info for link with subroutine sacio
c
c   
c scaio is in Tom Owens' SAC input/output routines in Subs library'
c	
      real instr
      integer year,jday,hour,min,isec,msec,itype
      character*8 sta,cmpnm,evnm
      character  idcode*2,staname*4,clenth*8,cpname*1
      common /tjocm/ dmin,dmax,dmean,year,jday,hour,min,isec,msec,sta,
     *            cmpnm,az,cinc,evnm,baz,delta,rayp,depth,decon,agauss,
     *              c,tq,instr,dlen,begin,t0,t1,t2
c
c **************************************************************************
c
c   parameter definitions may be found in sacio comments
c
      common /win/ data
      common /innout/ inunit,ounit
      inunit=5
      ounit=6
      zero1=cmplx(0.,0.)
      pi=3.141592654
      call iniocm
      call asktxt('Specify SAC file: ',eqfile)
      print *, eqfile
      call asktxt('Specify output file: ',outfil)
      print *, "itype input:(1-Sp 2-S2s1p 3-S1s2p 4-S3p)"
      read (*,*) itype
      print *, itype
      call sacio(eqfile,d2(1),npts,dt,1)
      write(*,*)delta,slowness
      call compslow(delta,slowness)
c      itype=1
      write(*,*)delta,slowness
      call psmoutt(d2,npts,begin,dt,slowness,itype)
      call sacio(outfil,d2(1),npts,dt,-1)
      end



      subroutine compslow(delta,slow)
      
        dimension x(36),y(36)
        data x/
     #  30.000,  32.000,  34.000,  36.000,  38.000,  40.000,  42.000,
     #  44.000,  46.000,  48.000,  50.000,  52.000,  54.000,
     #  56.000,  58.000,  60.000,  62.000,  64.000,  66.000,  68.000,
     #  70.000,  72.000,  74.000,  76.000,  78.000,  80.000,  82.000,
     #  84.000,  86.000,  88.000,  90.000,  92.000,  94.000,  96.000,
     #  98.000, 100.000/
        data y/
     #  15.664,  15.556,  15.431,  15.283,  15.121,  14.946,  14.760,
     #  14.567,  14.367,  14.160,  13.950,  13.734,  13.519,
     #  13.300,  13.079,  12.855,  12.630,  12.404,  12.176,  11.945,
     #  11.713,  11.477,  11.240,  11.000,  10.757,  10.509,  10.256,
     #   9.999,   9.734,   9.465,   9.186,   8.903,   8.699,   8.593,
     #   8.427,   8.323/


      if(delta.lt.x(1))then
      slow=y(1)-(y(2)-y(1))/(x(2)-x(1))*(x(1)-delta)
      end if
      if(delta.gt.x(36))then
      slow=y(36)+(y(36)-y(35))/(x(36)-x(35))*(delta-x(36))
      end if

      do 100 i=2,36
      if(delta.gt.x(i-1).and.delta.le.x(i))then
      slow=y(i-1)+(y(i)-y(i-1))/(x(i)-x(i-1))*(delta-x(i-1))
      end if
100   continue
      return
      end

       


C  psmout.f
C
C  The program makes moveout corrections to Ps delay times relative to P
C  onset for traces with different slownesses. It can also correct the
C  moveout for other 3 multiple phases optionally, they are Ppps, Ppss
C  and Psss.
C  Reference slowness is 6.4 s/deg.
C  It is assumed that Ps and the other multiples have the same slowness
C  with the corresponding P phase which may have big deviations for 
C  multiples at large depth.
C  It uses Earth-Flattening-Approximation.
C
C  Input data are ASCII written from SeismicHandler with the header: 
C  SLOWNESS" (TMP_PSMOUT.ASC)
C
C  Model: ./iasp91.dat
C
C  Output: binary file of SeismicHandler.
C  (PSMOUT.QBN, PPPSMOUT, PPSSMOUT, PSSSMOUT)
C
C  Xiaohui Yuan, Sep. 1997.
C
C
        subroutine psmoutt(x,npt,tbeg,delt,slowness,itype)
	dimension ttmo(100,200),tt(200),maxt(100),maxy(2000)
        dimension x(5000),y(5000),ysum(5000),tref(5000)
	character name*20,ch*10
	common /slow/ slow1,slow2,dslow
	common /time/ time1,time2,dtime

c
c  i/o channels
c
	in=31
	out=32
c
c  open input and output data file and read some values
c
c        open(in,file='TMP_PSMOUT.ASC',status='old')
c        read(in,'(a)') name
c        read(in,*) ntrace,tbeg,tend,delt,itype

        ntrace=1
        tend=tbeg+(npt-1)*delt
c        numt=(tend-tbeg)/delt+1
        numt=npt
	do 4 j=1,numt
4	tref(j)=tbeg+(j-1)*delt

c
c   open output ascii files for different phase type:
c     itype=1: Ps
c     itype=2: P2p1s
c     itype=3: P1p2s
c     itype=4: P3s
c
C23456789012345678901234567890123456789012345678901234567890123456789012
c        if (itype.eq.1) then
c           open(out,file='PSMOUT.QBN',status='unknown'
c     &             ,form='unformatted',access='direct',recl=4)
c        else if (itype.eq.2) then
c           open(out,file='PPPSMOUT.QBN',status='unknown'
c     &             ,form='unformatted',access='direct',recl=4)
c        else if (itype.eq.3) then
c           open(out,file='PPSSMOUT.QBN',status='unknown'
c     &             ,form='unformatted',access='direct',recl=4)
c        else if (itype.eq.4) then
c           open(out,file='PSSSMOUT.QBN',status='unknown'
c     &             ,form='unformatted',access='direct',recl=4)
c        else
c           write(*,*) 'error: wrong phase type'
c           stop
c        endif
c
c  moveout table
c
	call table(itype,nslow,maxt,ttmo)
c
	do 6 j=1,numt
6	ysum(j)=0
	do 8 i=1,ntrace
8	maxy(i)=numt
c
c  Loop for traces
c
	maxtr=0
	do 100 itr=1,ntrace
	write(*,96)itr
96	format('  trace: ',i5)

c
c  read in data of one trace
c
c        read(in,*) ch,delta
c        read(in,*) ch,length
c        read(in,*) ch,slowness
c        read(in,*) (x(j),j=1,length)
c        read(in,*)

        delta=delt
        length=npt
	if (slowness.gt.slow2) then
	   write(*,89) slow2
89	   format(' slowness exceeds',f10.3,', neglected')
	   goto 100
	endif
	maxtr=maxtr+1
c
c  moveout table for this trace
c
	is1=int((slowness-slow1)/dslow)+1
	if(is1.eq.nslow) is1=is1-1
	is2=is1+1
	si1=slow1+(is1-1)*dslow
	si2=si1+dslow
	maxtable=maxt(is2)
	do 10 j=1,maxtable
	tt(j)=ttmo(is1,j)+(ttmo(is2,j)-ttmo(is1,j))
     &        *(slowness-si1)/(si2-si1)

10	continue
c
c  moveout correction
c
	do 20 j=1,numt
	if (tref(j).le.0) then
	   treal=tref(j)
	else
	   i1=int((tref(j)-time1)/dtime)+1
	   i2=i1+1
	   if (i2.gt.maxtable) then
	      write(*,911) treal
911	      format('  warning: maximum Sdp delay time: ',f10.3,'s')
	      maxy(itr)=j-1
	      goto 21
	   endif
	   ti1=time1+(i1-1)*dtime
	   ti2=ti1+dtime
	   ttx=tt(i1)+(tref(j)-ti1)*(tt(i2)-tt(i1))/(ti2-ti1)
	   treal=tref(j)+ttx
	end if

	i1=int((treal-tbeg)/delta)+1

	if (i1.gt.length) then
	   write(*,*) '  warning: data not enough length'
	   maxy(itr)=j-1
	   goto 21
	else if (i1.eq.length) then
	   i1=i1-1
	else
	   i2=i1+1
	   ti1=tbeg+(i1-1)*delta
	   ti2=ti1+delta
	   y(j)=x(i1)+(treal-ti1)*(x(i2)-x(i1))/(ti2-ti1)
	   ysum(j)=ysum(j)+y(j)
	endif
20	continue
21	continue
c
c  output to file
c
	do 30 j=1,numt
	if (j.gt.maxy(itr)) y(j)=0
30	continue
c************************************
c        do 40 j=1,numt
c        recnum=j+(itr-1)*numt
c40      write(out,rec=recnum) y(j)
c************************************
        do 40 j=1,numt
40      x(j)=y(j)

100	continue

        return
	end


********************************************************************

	subroutine table(itype,nslow,maxt,ttpsmo)
	common /slow/ slow1,slow2,dslow
	common /time/ time1,time2,dtime
	dimension depco(500),vpco(500),vsco(500),layco(500)
	dimension dep(2000),vp(2000),vs(2000),hk(2000),ro(2000)
	dimension slowness(100),depth(200),ttps(100,200),ttps0(200)
	dimension t0(200),ttpsmo(100,200),maxdep(100),maxt(100)
	character model*50
c
c  parameters for moveout table
c
	slow1=0
	slow2=14
	dslow=0.2
	nslow=int((slow2-slow1)/dslow)+1

	depth1=0
	depth2=2000
	ddepth=20
	ndepth=int((depth2-depth1)/ddepth)+1

	time1=0
	time2=200
	dtime=2
	ntime=int((time2-time1)/dtime)+1

	do 3 i=1,nslow
3	maxt(i)=ntime

c
c  read in velocity model
c
	model='iasp91.dat'
	open(11,file=model,status='old')
	read(11,*)
	read(11,*)
	nlayco=0
	do 10 i=1,500
	read(11,*) depco(i),vpco(i),vsco(i),layco(i)
	if (abs(depco(i)).lt.0.001.and.abs(vpco(i)).lt.0.001
     1	.and.abs(vsco(i)).lt.0.001.and.abs(layco(i)).lt.0.001) 
     2	goto 11
	nlayco=nlayco+1
10	continue
11	close(11)

	maxdep0=0
	do 301 i=1,nslow
301	maxdep(i)=0

c
c  Loop for depth
c
	do 200 ldep=1,ndepth
	depthmax=depth1+ddepth*(ldep-1)
	depth(ldep)=depthmax
c
c  transformation of gradient model to homogeneous layered model
c
	lay=0
	do 20 i=1,nlayco
	if(layco(i).eq.0) goto 20
	vp0=vpco(i-1)
	vs0=vsco(i-1)
	dep0=depco(i-1)
	devp=(vpco(i)-vpco(i-1))/layco(i)
	devs=(vsco(i)-vsco(i-1))/layco(i)
	dedep=(depco(i)-depco(i-1))/layco(i)
	do 19 j=1,layco(i)
	lay=lay+1
	vp(lay)=vp0+devp*(j-0.5)
	vs(lay)=vs0+devs*(j-0.5)
	hk(lay)=dedep
	dep(lay)=dep0+dedep*(j-1)
	depm=dep(lay)+dedep
	if (depm.ge.depthmax) then
	   hk(lay)=hk(lay)-(depm-depthmax)
	   goto 21
	end if
19	continue
20	continue
21	continue
c
c  Earth Flattening Approximation
c
	call efa(lay,dep,hk,vp,vs,ro,1)
c
c  for reference slowness p=6.4
c
	slow=6.4/111.195
	ttp=0
	tts=0
	xp=0
	xs=0
	do 80 i=1,lay
	pvp=slow*vp(i)
	pvs=slow*vs(i)

c   judge if post critical angle
	if (pvp.ge.1) goto 201

	h=hk(i)
	termp=1/sqrt(1-pvp*pvp)
	terms=1/sqrt(1-pvs*pvs)
	ttp=ttp+h*termp/vp(i)
	tts=tts+h*terms/vs(i)
	xp=xp+h*pvp*termp
	xs=xs+h*pvs*terms
80	continue
c
c   computes delay time for different phase type:
c     itype=1: Ps
c     itype=2: P2p1s
c     itype=3: P1p2s
c     itype=4: P3s
c
	if (itype.eq.1) then
	   ttps0(ldep)=tts-ttp-(xs-xp)*slow
	else if (itype.eq.2) then
	   ttps0(ldep)=tts+ttp-(xs+xp)*slow
	else if (itype.eq.3) then
	   ttps0(ldep)=2*tts-2*xs*slow
	else if (itype.eq.4) then
	   ttps0(ldep)=3*tts-ttp-(3*xs-xp)*slow
	else
	   write(*,*) 'error in subroutine: wrong phase type'
	   stop
	endif

	maxdep0=maxdep0+1

c
c  Main calculation, loop for slowness
c
	do 101 lslo=1,nslow
	slowness(lslo)=slow1+dslow*(lslo-1)
	slow=slowness(lslo)/111.195
	ttp=0
	tts=0
	xp=0
	xs=0

	do 100 i=1,lay
	pvp=slow*vp(i)
	pvs=slow*vs(i)

c   judge if post critical angle
	if (pvp.ge.1) goto 101

	h=hk(i)
	termp=1/sqrt(1-pvp*pvp)
	terms=1/sqrt(1-pvs*pvs)
	ttp=ttp+h*termp/vp(i)
	tts=tts+h*terms/vs(i)
	xp=xp+h*pvp*termp
	xs=xs+h*pvs*terms
100	continue
c
c   computes delay time for different phase type:
c     itype=1: Ps
c     itype=2: P2p1s
c     itype=3: P1p2s
c     itype=4: P3s
c
	if (itype.eq.1) then
	   ttps(lslo,ldep)=tts-ttp-(xs-xp)*slow
	else if (itype.eq.2) then
	   ttps(lslo,ldep)=tts+ttp-(xs+xp)*slow
	else if (itype.eq.3) then
	   ttps(lslo,ldep)=2*tts-2*xs*slow
	else if (itype.eq.4) then
	   ttps(lslo,ldep)=3*tts-ttp-(3*xs-xp)*slow
	else
	   write(*,*) 'error in subroutine: wrong phase type'
	   stop
	endif
	maxdep(lslo)=maxdep(lslo)+1
101	continue
200	continue
201	continue
c
c  differences of delay times
c
	do 120 i=1,nslow
	do 120 j=1,maxdep(i)
	ttps(i,j)=ttps(i,j)-ttps0(j)
120	continue
c
c  Interpolation for an equal-sampling reference delay-time array.
c
	do 130 i=1,nslow

	ind=1

	do 129 j=1,ntime
	t0(j)=time1+(j-1)*dtime

141	continue
	if (ttps0(ind).ge.t0(j)) then
	   if (ind.eq.1) then
	      ttpsmo(i,j)=ttps(i,1)
	   else if (ind.gt.maxdep(i)) then
	      maxt(i)=j-1
	      goto 130
	   else
	      ind1=ind-1
	      ttpsmo(i,j)=ttps(i,ind1)+(t0(j)-ttps0(ind1))*(ttps(i,ind)
     &                    -ttps(i,ind1))/(ttps0(ind)-ttps0(ind1))
	   end if
	else
	   ind=ind+1
	   goto 141
	end if

129	continue
130	continue

	return
	end

************************************************************

      subroutine efa(n,dep,h,vp,vs,ro,i)
c
c--- earth flattenning approximation for gradient model
c--- i=1(from spherical to flat), i.ne.1 (inverse)
c--- modefied from kosarev's "efal" in "bip.f"
c
      dimension dep(n),h(n),vp(n),vs(n),ro(n)
      if(i.eq.1) then
         za=0.
         zz=0.
         do k=1,n
c      write(*,*)'k=',k,' y1',y1,' dep,vp,vs,ro,h'
c     &,dep(k),vp(k),vs(k),ro(k),h(k)
              z1k=zz+h(k)*0.5
              zz=zz+h(k)
              y1=6370./(6370.-z1k)
              yz=6370./(6370.-zz)
              zm=6370.*alog(yz)
              h(k)=zm-za
              za=zm
              vp(k)=vp(k)*y1
              vs(k)=vs(k)*y1
              ro(k)=ro(k)/y1
              if (k.gt.1) dep(k)=dep(k-1)+h(k-1)
c      write(*,*)'k=',k,' y1',y1,' dep,vp,vs,ro,h'
c     &,dep(k),vp(k),vs(k),ro(k),h(k)
         end do
      else
         za=0.
         zz=0.
         do k=1,n
              z1k=zz+h(k)*0.5
              if(k.eq.n) z1k=zz
              zz=zz+h(k)
              z1k=6370.*(1.-exp(-z1k/6370.))
              zm=6370.*(1.-exp(-zz/6370.))
              y1=6370./(6370.-z1k)
              h(k)=zm-za
              za=zm
              vp(k)=vp(k)/y1
              vs(k)=vs(k)/y1
              ro(k)=ro(k)*y1
              if (k.gt.1) dep(k)=dep(k-1)+h(k-1)
         end do  
       end if
       return
       end
