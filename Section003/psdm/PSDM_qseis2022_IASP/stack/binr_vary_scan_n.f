C main program for Midpoint Binning and Moveout Corrections
C main program
c  rectangle bin!!!!!!!!!!!!!!!!!!!!!!
C* revised by Ling Chen for varying bin size, 10/30/2003
C* revised again for multi var ranges for stacking  11/26/2003
c  no time-->depth
C*
C* revised from binr_vary_new_r.f on 08/05/04:
C* minimum YBIN, DYBIN, maximum YBIN and XBIN are all depth-dependent
C*
C* revised from binr_vary_new_rr.f on 08/08/04:
C* multiple-profile stacking by scanning the end-point locations and trending
C* revised from binr_vary_scan.f on 01/21/05:
C* reading multiple piercing point data files and stacking
C* 
C* Modified further by calling subroutine utm_geo to calculate
C* the geographic locations of the profiles
C* Modified to read UTM_PROJECTION_ZONE from the input file by Ling
C* -- 01/29/14, IGGCAS
*$ debug
	 include "stack.inc"
	 parameter(ism0=5,ism1=2)
	 common/param1/baz(num),gcarc(num),stla(num),stlo(num)
       common/param2/ppo(num),icnt(num),n,limg,npief0
	 common/param3/piercx(num,nwmax),piercy(num,nwmax),tto(num)
       COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC),NW,NWI(nwmax),
     & NWID(nwmax)
	 common/locate1/begla,beglo,endla,endlo,stepbin,SA,CA,leasttra,
     & XBIN(nwmax),YBIN0(nwmax),DYBIN(nwmax),YBINM(nwmax)
	 common/locate2/numsta,numtra,numbin,dep(nwmax),
     & avergca,averaz	 
       common/param4/nsta(ina,nwmax),ntra(ina,nwmax),depth(icc),
     & YBIN(ina,nwmax)
	 common/param5/rang(icount),pp(icount)
	 common/param6/pdstime(icount,icc),pdst0(icc)
	 common/param7/dpdstime(icount,icc),dpdst0(icc)
       dimension x(ina),y(ina),stack(0:maxnump),temp5(0:maxnump)
       dimension idepb(nwmax),idepe(nwmax),stack1(0:maxnump)
	 character*150 iaj,ibj,pierc_data,inname,inn1,inn2
	 character*150 inname1(num),inname2(num),inname3
	 character*150 DIR,DIRimg,outfile,timefile
	 integer inw0(nwmax),npt,np0,idist
	 integer lxbin,lybin0,ldybin,lybinm,lnumtra,linw0,linw1,ldx,ly
	 character*10 cxbin,cybin0,cdybin,cybinm,cnumtra,cinw0,cinw1,
     & cdx,cy,cnorm
	 integer lname(5),lbiname(maxnum)
	 real t00,t01
	 integer ninw,inw,inw1,it0,it1,itt,itb,imethod,inorm,lnorm
	 integer nsp,nalp,noutd,ioutd(nwmax),istack,lenout,lensp
	 real xlenp,XBIN0(nwmax)
	 character*20 csp
	 integer npief,ioutb                         ! 01/21/05
	 double precision wx0,wy0,wx,wy,wlon,wlat    ! new
	 integer UTM_PROJECTION_ZONE,iway            ! new
	 logical SUPPRESS_UTM_PROJECTION             ! new

*CL new: parameters used in subroutine utm_geo
	 SUPPRESS_UTM_PROJECTION = .false.
*01/29/14	 UTM_PROJECTION_ZONE = 50                    ! this number may be different for different regions

c---------------
	 cont=pi/180.
	 open(7,file='binr_vary_scan_n.inp',status='old')
	 open(8,file='bin.out',status='unknown')
*<1>
	 write(*,*)'pl. read begin,end coordinate of starting points (km),
     & and point interval (km)'
	 read(7,'(a)') ibj
	 read(7,*) begla0,beglo0,endla0,endlo0,dsp
	 write(*,*) begla0,beglo0,endla0,endlo0,dsp
	 write(8,*) 'begin and end coordinate of starting points (km),
     & and point interval (km)'
	 write(8,*) begla0,beglo0,endla0,endlo0,dsp
*<1>
	 write(*,*)'pl. read length of the profiles(km) and azimuth range'
	 read(7,'(a)') ibj
	 read(7,*) xlenp,alphab,alphae,dalp
	 write(*,*) xlenp,alphab,alphae,dalp
	 write(8,*)'length of the profiles(km) and azimuth range (degree)'
	 write(8,*) xlenp,alphab,alphae,dalp
	 alphab=alphab*pi/180.0
	 alphae=alphae*pi/180.0
	 dalp=dalp*pi/180.0
*<2>
	 write(*,*) 'pl. read the spacing between bins (km),
     & and least number of traces, ratio of least number: rnumtra
     & and UTM_PROJECTION_ZONE'
	 read(7,'(a)') ibj
	 read(7,*) stepbin,numtra,rnumtra,UTM_PROJECTION_ZONE !01/29/14               ! 08/06/04
	 if(rnumtra.le.0.)rnumtra=1.0                   ! 03/18/05
	 write(*,*) stepbin,numtra,rnumtra              ! 08/06/04
	 write(8,*) 'the spacing between bins (km),least number of traces'
	 write(8,*) stepbin,numtra,rnumtra
	 numtra0=numtra*rnumtra
*<3>
	 write(*,*)'pl. read the time file name for stacking: timefile'
	 read(7,'(a)') ibj
	 read(7,'(a)') timefile
	 write(8,*)'the time file name: timefile'
	 write(8,'(a)') timefile
*<4>
	 write(*,*)'pl. read the output file name for stacking: outfile'
	 read(7,'(a)') ibj
	 read(7,'(a)') outfile
	 write(8,*)'the output file name: outfile (no postfix)'
	 write(8,'(a)') outfile
*<5>
	 write(*,*) 'pl. read ouput number of time samples in each trace'
	 read(7,'(a)') ibj
	 read(7,*) npt,dt
	 write(*,*) npt,dt
	 write(8,*) 'ouput number of time samples in each trace, dt'
	 write(8,*) npt,dt
*<6>
	 write(*,*) 'pl. read the indexes of reference ray among 1 -- nw
     & and corresponding XBIN'
	 read(7,'(a)') ibj
CL: 12/26/03: consider multiple focusing depths
	 read(7,*) ninw, (inw0(i),XBIN(i),i=1,ninw)
* XBIN(1) < 0: same XBIN for all depths
	 if(XBIN(1).lt.0.)then
	    XBIN(1)=-1.*XBIN(1)
	    do i=2,ninw
	      XBIN(i)=XBIN(1)
	    enddo
	 endif
	 do i=1,ninw
	   XBIN0(i)=XBIN(i)
	 enddo
C* Warning: here may need to sort the elements of inw in asccending way
	 write(*,*) 'inw =',(inw0(i),XBIN(i),i=1,ninw)
	 write(*,*)
	 write(8,*) 'the indexes of reference ray among 1 -- nw'
	 write(8,*) ninw,(inw0(i),XBIN(i),i=1,ninw)
*<7>
	 write(*,*) 'pl. read the corresponding minimum YBIN'
	 read(7,'(a)') ibj
	 read(7,*) (YBIN0(i),i=1,ninw)
* YBIN0(1) < 0: same YBIN0 for all depths
	 if(YBIN0(1).lt.0.)then
	    YBIN0(1)=-1.*YBIN0(1)
	    do i=2,ninw
	      YBIN0(i)=YBIN0(1)
	    enddo
	 endif
	 write(*,*) (YBIN0(i),i=1,ninw)
	 write(*,*)
	 write(8,*) 'the corresponding minimum YBIN'
	 write(8,*) (YBIN0(i),i=1,ninw)
*<7>
	 write(*,*) 'pl. read the corresponding DYBIN'
	 read(7,'(a)') ibj
	 read(7,*) (DYBIN(i),i=1,ninw)
* DYBIN(1) < 0: same DYBIN for all depths
	 if(DYBIN(1).lt.0.)then
	    DYBIN(1)=-1.*DYBIN(1)
	    do i=2,ninw
	      DYBIN(i)=DYBIN(1)
	    enddo
	 endif
	 write(*,*) (DYBIN(i),i=1,ninw)
	 write(*,*)
	 write(8,*) 'the corresponding DYBIN'
	 write(8,*) (DYBIN(i),i=1,ninw)
*<7>
	 write(*,*) 'pl. read the corresponding maximum YBIN'
	 read(7,'(a)') ibj
	 read(7,*) (YBINM(i),i=1,ninw)
* YBINM(1) < 0: same YBINM for all depths
	 if(YBINM(1).lt.0.)then
	    YBINM(1)=-1.*YBINM(1)
	    do i=2,ninw
	      YBINM(i)=YBINM(1)
	    enddo
	 endif
	 write(*,*) (YBINM(i),i=1,ninw)
	 write(*,*)
	 write(8,*) 'the corresponding maximum YBIN'
	 write(8,*) (YBINM(i),i=1,ninw)
*<8>
	 write(*,*) 'pl. read temporary directory name to store the
     & intermedial files (.img)'
	 read(7,'(a)') ibj
	 read(7,'(a)') DIRimg
	 write(*,'(a)') DIRimg
*<9>
	 write(*,*)'pl. read the moveout index idist'
* idist= 0: using gcarc0 in pierc_new.dat for moveout; <0: no moveout; 
*        1: selected distance; else: time shift to p=0 approximately)'
	 read(7,'(a)')ibj
	 read(7,*)idist,gcarc1
* gcarc1 is used only when idist=1
	 write(*,*) idist,gcarc1
	 write(8,*)'the moveout index idist =',idist
	 if(idist.eq.1)write(8,*)'reference gcarc =',gcarc1
*<10>
	 write(*,*) 'pl. read the mornalization index inorm'
* inorm = 0: no normalization; else: normalization for individual RFs
	 read(7,'(a)') ibj
	 read(7,*) inorm
	 write(*,*) 'inorm =',inorm
	 write(8,*) 'normalization index: inorm'
	 write(8,*) inorm
*<11>
	 write(*,*) 'pl. read the output number and depth indexes in ninw'
	 read(7,'(a)') ibj
	 read(7,*) noutd,(ioutd(i),i=1,noutd)
	 write(*,*) noutd,(ioutd(i),i=1,noutd)
	 write(8,*) 'output number and depth indexes in ninw'
	 write(8,*) noutd,(ioutd(i),i=1,noutd)
*<12>
	 write(*,*) 'pl. read the output index for stacking: istack'
* istack = 0: not perform stacking;
*       else: perform stacking and output stackinf results
	 read(7,'(a)') ibj
	 read(7,*) istack
	 write(*,*) 'stacking index: istack =',istack
	 write(8,*) 'index for stacking: istack =',istack
* end of parameter input

c-------------------
* 01/21/05:
*<13>
	 write(*,*) 'pl. read the output index for baz: ioutb'
* ioutb = 0: not output back-azimuth;
*      else: output back-azimuth
	 read(7,'(a)') ibj
	 read(7,*) ioutb
	 write(*,*) 'output index for baz: ioutb =',ioutb
*<14>
	 write(*,*)'piercing point data file number npief =?'
	 read(7,'(a)') ibj
	 read(7,*) npief0                                     ! 04/21/06
	 npief=iabs(npief0)
	 write(*,*) 'npief =',npief
	 write(8,*) 'piercing point data file number npief =',npief
	 write(*,*)'piercing point data file name =?'
	 read(7,'(a)') ibj
	 i0=0
	 gca=0.0
	 do 1 ipe=1,npief
	   read(7,'(a)') pierc_data
	   write(*,'(a)') pierc_data

c	 read pierc_data
	   open(9,file=pierc_data,status='unknown')
C* nw and NWI should be the same in different piercing point data files,
C* and same for igca
	   read(9,*) nw,(NWI(I1),dep(i1),i1=1,nw)
	   read(9,*) evla0,evlo0
C* new:
	   wlon=dble(evlo0)
	   wlat=dble(evla0)
	   iway = 0
	   call utm_geo(wlon,wlat,wx0,wy0,UTM_PROJECTION_ZONE,iway,
     ^SUPPRESS_UTM_PROJECTION)
C* new end
	   read(9,*) gcarc0,pp0,aminp,averaz
	   read(9,*) n,np0
* np0 is the time sampling point in pierc_new.dat, no meaning
	   read(9,'(a)')DIR
	   ldir=lengths(DIR)
	   do 10 i=i0+1,i0+n
         read(9,'(a)') inname
         read(9,*) icnt(i),stla(i),stlo(i)
111	 format(1x,a30,i5,2f15.0)
	   linname=lengths(inname)
	   inname1(i)=DIR(1:ldir)//inname(1:linname)
* 12/16/03: considering multi-direction cases
	   lndir=0
	   do j=linname,1,-1
	     if(inname(j:j).eq.'/') then
             lndir=lndir+1
             lname(lndir)=j
	     endif
	   enddo
	   do j=1,lndir-1
	     j1=lname(j)
	     inname(j1:j1)='_'
	   enddo
	   j1=lname(lndir)
109	   inname2(i)=DIRimg(1:limg)//inname(j1+1:linname-3)//'img'
	   read(9,112) gcarc(i),baz(i),cosazb,ppo(i)
112	 format(1x,2f8.2,2e12.4)
	   read(9,113) (piercx(i,i2),i2=1,nw),(piercy(i,i2),i2=1,nw)
113	 format(1x,10f8.2)
	   read(9,114) (tto(i2),i2=1,nw)
114	 format(1x,10f8.2)
10	   continue	  	  	  
         close(9)

	   i0=i0+n
	   gca=gca+gcarc0*float(n)
1	 continue
	 n=i0
	 gcarc0=gca/float(n)
	 write(*,*)'RF number =',n,' gcarc0 =',gcarc0

	 close(7)
* end 01/21/05
	
	limg=lengths(DIRimg)
	inname3=DIRimg(1:limg)//'/test.img'
	
c--------------------
c	 claculating the number of starting points, azimuths and stack bins
	 nsp=int(sqrt((endla0-begla0)**2+(endlo0-beglo0)**2)/dsp+0.001)+1
	 write(8,*) 'the number of starting points:  nsp =',nsp
	 write(*,*) 'the number of starting points:  nsp =',nsp
	 write(*,*)begla0,endla0,beglo0,endlo0,dsp
	 if(abs((endlo0-beglo0)).lt.0.1) then
	 alphap=pi/2.
	 else
	 alphap=atan((endla0-begla0)/(endlo0-beglo0))
	 endif
	 write(8,*) 'the alphap is:',alphap
	 SAP=SIN(alphap)
	 CAP=COS(alphap)
	 WRITE(*,*)  'the alphap is',ALPHAP,SAP,CAP
	 nalp=int((alphae-alphab)/dalp+0.001)+1
	 write(8,*) 'the number of azimuth: nalp =',nalp
	 write(*,*) 'nalp =',nalp
	 numbin=int(xlenp/stepbin)+1
	 write(8,*) 'the number of stack bin: numbin =',numbin
	 write(*,*) 'numbin =',numbin

	lenout=lengths(outfile)
	inn1=outfile(1:lenout)//'_num.dat'
	nrecl=4*numbin*noutd
      open(20,file=inn1,status='unknown',recl=nrecl,
     &	access='direct')
	inn2=outfile(1:lenout)//'_yb.dat'
      open(21,file=inn2,status='unknown',recl=nrecl,
     &	access='direct')
	inname=outfile(1:lenout)//'_profile.txt'
	open(22,file=inname,status='unknown')
	write(8,*) 'output file name containing RF numbers in the bins:'
	write(8,'(a)')inn1
	write(8,*) 'output file name containing YBIN of the bins:'
	write(8,'(a)')inn2
	write(8,*) 'output file name containing start and end points:'
	write(8,'(a)')inname

C* 01/21/05
	if(ioutb.ne.0) then
* for dep = 32 (cwbq), gcarc = 30 => r = 9.5 km, gcarc = 90 => r = 4.75
	  slope=(9.5-4.75)/(90.-30.)
	  gca1b=4.75+slope*(0-30)          ! 2.375
	  gca1e=4.75+slope*(100-30)        ! 10.2917
	  rat=10./(gca1e-gca1b) 
	  nrecl=4*6
	  inn2=outfile(1:lenout)//'_baz.dat'
        open(30,file=inn2,status='unknown',recl=nrecl,
     &	access='direct')
	  do i=1,n
	    gca1=(4.75+slope*(gcarc(i)-30.)-gca1b)*rat
	    gca_bazx=gca1*cos(baz(i)*cont)
	    gca_bazy=gca1*sin(baz(i)*cont)
	    write(30,rec=i)gcarc(i),baz(i),piercx(i,1)-stla(i),
     & piercy(i,1)-stlo(i),gca_bazx,gca_bazy
	  enddo
	  close(30)
	endif

	rlat=1./rad

	if(istack.ne.0) then
c-------------------
c	read moveout parameter 
	 open(11,file=timefile,status='unknown')
	 read(11,*) nre,nrange
	 read(11,*) (depth(i1),i1=1,nre)	
	 do 30 i=1,nrange
	 read(11,*) rang(i),pp(i)
	 read(11,*) (pdstime(i,j),j=1,nre)
30	continue	 
	 close(11) 
c	 write(*,*) 'nre=',nre,'   nrange=',nrange
c	 write(*,*) depth(10),rang(10),pp(10),pdstime(10,10)

c	 find the average rang(i) and corresponding pdst0(j)
	 if(idist.eq.1)gcarc0=gcarc1
	 ddmin=10.
	 do 40 i=1,nrange
	 ddmin1=abs(gcarc0-rang(i))
	 if(ddmin1.lt.ddmin) then
	 ddmin=ddmin1
	 i0=i
	 endif
40	 continue
         avergca=rang(i0)      ! get the average epi. distance
         do 50 i=1,nre
         pdst0(i)=pdstime(i0,i)
	 dpdst0(i)=pdstime(i0,i)-pdstime(nrange,i)
50       continue
         write(*,*) avergca

c--------------------
c	 do moveout correction and time-depth migration
	   write(*,*) 'caution! all record  must have the same delt'
	   nrecl=4*npt
	   imgout=13
         open(imgout,file=inname3,status='unknown',recl=nrecl,
     &	access='direct')
	  do 60 i=1,n
	   ip=i
	   call moveout(inname1,imgout,ip,nre,nrange,npt,idist,inorm)
60	  continue
	  close(imgout)

	  nrecl=4*npt
        open(17,file=inname3,status='old',recl=nrecl,
     &	access='direct')

	itt=-1
	do 495 inwd=1,ninw
	   inw=inw0(inwd)
	   t00=pdst0(nwi(inw))
	   if(inwd.lt.ninw)then
	     inw1=nwi(inw0(inwd+1))
	   else
	     inw1=nre
	   endif
	   t01=pdst0(inw1)
	   if(idist.gt.1) then
	     t00=t00-dpdst0(nwi(inw))
	     t01=t01-dpdst0(inw1)
	   endif
	   it0=int(t00/dt+0.5)
	   it1=int(t01/dt+0.5)
	   idepb(inwd)=itt+1
	   itt=it0+(it1-it0)*0.4
	   itt=min(npt-1,itt)
	   if(inwd.eq.ninw)itt=npt-1
	   idepe(inwd)=itt
	   write(*,*)nwi(inw),t00,inw1,t01
	   write(*,*)'inw =',inw,' itb =',idepb(inwd),' ite =',idepe(inwd)
495	continue

	   write(8,*)'output stacked RF data file:'
	endif          ! endif(istack.ne.0)
C****************************************************************************

	 do 1000 isp=1,nsp
       write(*,*)'isp =',isp

	 begla=begla0+real(isp-1)*dsp*sin(alphap)
       beglo=beglo0+real(isp-1)*dsp*cos(alphap)
	 ispp=(isp-1)*nalp
	 
	 if(istack.ne.0) then
	   nrecl=4*npt
	   if(nsp.gt.1) then
	     call num2str(isp,csp,lensp)
	     inname='stack_'//outfile(1:lenout)//'_'//csp(1:lensp)//'.dat'
	   else
	     inname='stack_'//outfile(1:lenout)//'.dat'
	   endif
         open(15,file=inname,status='unknown',recl=nrecl,
     &	access='direct')
	   write(8,'(a)') inname
	 endif

	 do 999 ialp=1,nalp
	 irecd=ispp+ialp
	 irecd1=(ialp-1)*numbin

	 alpha_r=alphab+real(ialp-1)*dalp
	 alpha=pi*0.5-alpha_r
	 SA=SIN(alpha)
	 ca=cos(alpha)
*	 WRITE(*,*)  'the alpha is',ALPHA,SA,CA

        do i=1,numbin
	    inumb=i
           x1=begla+real(i-1)*stepbin*sin(alpha)
           y1=beglo+real(i-1)*stepbin*cos(alpha)
           x(inumb)=x1
           y(inumb)=y1
*	    write(8,*) 'the number of bin=',inumb,x1,y1
	 enddo
	 if(nsp*nalp.gt.10)then
	   write(22,*)y(1),x(1),y(numbin),x(numbin)
	 else
C* new
	   iway = 1
	   do i=1,numbin
	      wx=wx0+dble(y(i)*1000.)
	      wy=wy0+dble(x(i)*1000.)
	      call utm_geo(wlon,wlat,wx,wy,UTM_PROJECTION_ZONE,iway,
     ^SUPPRESS_UTM_PROJECTION)
            write(22,1502)wlon,wlat,y(i),x(i),(i-1)*stepbin
*            write(22,1502)wlon,wlat,y(i),x(i)
*            write(*,1502)wlon,wlat,wx,wy
	   enddo
1502	format(f12.4,f12.5,2X,f15.4,f15.4,f12.1)
C* new end
	   write(22,*)
	 endif

c--------------------
CL: 12/26/03: consider multiple focusing depths
	nrecl=4*maxnum
      open(10,file='bin.dat',status='unknown',recl=nrecl,
     & access='direct')
         do 20 inumb=1,numbin
         call findrec(x(inumb),y(inumb),inumb,ninw,inw0)    ! 06/20/04
20	  continue
510	format(3x,i5,3X,I5,2X,f8.2,3X,I5,2X,f8.2)
      close(10)

	write(20,rec=irecd)
     & ((float(ntra(i,ioutd(j))),i=1,numbin),j=1,noutd)
	write(21,rec=irecd)((YBIN(i,ioutd(j)),i=1,numbin),j=1,noutd)

	if(istack.ne.0) then
c---------------------------
c        stack each bin to stack.dat
c-------------------

	   nrecl=4*maxnum
         open(16,file='bin.dat',status='unknown',recl=nrecl,
     &	access='direct')
         do 310 i=1,numbin
         do 320 j=0,npt-1
         stack(j)=0.
320      continue
	   i0=(i-1)*ninw
         do 330 j1=1,ninw
         read(16,rec=i0+j1)(lbiname(j3),j3=1,maxnum)

         do 340 j4=1,ntra(i,j1)
	     k=lbiname(j4)
1017     format(5f15.0)
         read(17,rec=k) (temp5(ji),ji=0,npt-1)
         do 350 j3=idepb(j1),idepe(j1)
         stack(j3)=stack(j3)+temp5(j3)
350      continue
340      continue

         if(ntra(i,j1).gt.numtra0) then                    ! 03/17/05
          do 370 j4=idepb(j1),idepe(j1)
           stack(j4)=stack(j4)/real(ntra(i,j1))
370       continue
         else
          do 380 j4=idepb(j1),idepe(j1)
           stack(j4)=stack(j4)/real(numtra0)               ! 03/17/05
380       continue
         endif

	 if(j1.gt.1)then
         it0=idepb(j1)-ism0
	   it0=max(0,it0)
	   it1=idepb(j1)+ism0
	   it1=min(it1,idepe(j1))
	   do 390 j4=it0,it1
	     j4b=j4-ism1
	     j4b=max(0,j4b)
	     j4e=j4+ism1
	     j4e=min(j4e,itt)
	     stack1(j4)=0.
	     nj4=j4e-j4b+1
	     do j5=j4b,j4e
	       stack1(j4)=stack1(j4)+stack(j5)
	     enddo
	     stack1(j4)=stack1(j4)/nj4
390	   continue
	   do 400 j4=it0,it1
	     stack(j4)=stack1(j4)
400	   continue
	 endif
330	 continue
	 write(15,rec=i+irecd1) (stack(j4),j4=0,npt-1)
310    continue
	 close(16)

	endif                 ! endif(istack.ne.0)

999	 continue
	 close(15)
1000	 continue
	close(20)
	close(21)
	close(22)
	
	write(*,*)'output data size in *num and *yb.dat:'
	write(*,*)numbin,' *',noutd,' *',nalp,' *',nsp
	write(8,*)'output data size in *num and *yb.dat:'
	write(8,*)numbin,' *',noutd,' *',nalp,' *',nsp
	if(istack.ne.0) then
	  close(17)
	  write(*,*)'output data size in each stacked RF file:',
     & npt,' *',numbin,' *',nalp
	  write(8,*)'output data size in each stacked RF file:',
     & npt,' *',numbin,' *',nalp
	endif
	close(8)

	stop
	end


	 
