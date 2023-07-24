	 subroutine moveout(inname1,imgout,ip,nre,nrange,np0,idist,inorm)
	 include "stack.inc"
	 parameter(nhead=158)
	 common/param1/baz(num),gcarc(num),stla(num),stlo(num)
       common/param2/ppo(num),icnt(num),n,ldir,npief0
	 common/param3/piercx(num,nwmax),piercy(num,nwmax),tto(num)
      COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC),NW,NWI(nwmax),
     &NWID(nwmax)
C	 common/locate1/begla,beglo,endla,endlo,stepbin,radii
	 common/locate1/begla,beglo,endla,endlo,stepbin,SA,CA,leasttra,
     & radiim(nwmax),dradi,radii0
	 common/locate2/numsta,numtra,numbin,dep(nwmax),
     & avergca,averaz	 
       common/param4/nsta(ina,nwmax),ntra(ina,nwmax),depth(icc),
     & radi(ina,nwmax)
	 common/param5/range(icount),pp(icount)
	 common/param6/pdstime(icount,icc),pdst0(icc)
	 common/param7/dpdstime(icount,icc),dpdst0(icc)
	 dimension temp1(maxnump),temp2(maxnump),temp3(maxnump)
	 dimension temp4(maxnump)
	 character*150 inname,inn1,inn2
	 character*150 inname1(num)
	 real tempp(nhead)
	 real dt(nwmax+1),amp0,amp
	 integer ndt(nwmax+1),mt(nwmax+1),move(nwmax+1),idist
	 real am_cor,pa
	
	lengthr=5
	nrecl=4*lengthr
C* 04/21/06
	if(npief0.lt.0) then                                      !wby 12/27/2010
      open(12,file=inname1(ip),status='old',recl=nrecl,
     &	access='direct',convert='BIG_ENDIAN')
	else
      open(12,file=inname1(ip),status='old',recl=nrecl,
     &	access='direct')
	endif
	 dista=10.
	 do 70 i3=1,nrange
	   dis1=abs(gcarc(ip)-range(i3))
         if(dis1.lt.dista) then
	     dista=dis1
	     id=i3
	   endif
70	 continue
* 09/19/04: amplitude correction
	pa=ppo(ip)
	pb=pa*pa
	am_cor = 151.5478*pb + 3.2896*sqrt(pb) + 0.2618 ! /* am corr for 3.0/3.5*/

 	 nt=int((pdstime(id,nre)-pdst0(nre))/0.1)
	 read(12,rec=1) delta,depmin,depmax,scale,odelta
	 read(12,rec=2) beg1,end1,omarker,amarker,ternal1
       npoint=nint((end1-beg1)*10./(delta*10.))+1
       nbeg=nint((0.-beg1)/delta)
       if(ip.eq.1) write(*,*)'total num =',npoint,' nbeg =',nbeg
115	 format(5f15.0)
116	 format(a80)
	 itemp=nhead/lengthr+1
	 ntotal=nhead+npoint
	 itotal=ntotal/lengthr
	 if(lengthr*itotal.lt.ntotal)itotal=itotal+1
	 itemp1=nhead-(itemp-1)*lengthr
	 itemp2=lengthr-itemp1
	 read(12,rec=itemp)(tempp(j1),j1=1,itemp1),(temp1(j1),j1=1,itemp2)
	 itemp2=itemp2+1
	 do j1=itemp+1,itotal
	   jj1=itemp2
	   jj2=jj1+4
	   jj2=min(jj2,npoint)
	   read(12,rec=j1)(temp1(j2),j2=jj1,jj2)
	   itemp2=itemp2+5
	 enddo
	 
	if(idist.ge.0) then
         do 85 i6=nbeg+1,npoint
           temp2(i6-nbeg)=temp1(i6)
85       continue
	   leng1=npoint-nbeg
         if(abs(nt).le.2) then
           do 90 j1=1,leng1
           temp3(j1)=temp2(j1)
90         continue
	     goto 180
         else
           if(nt.gt.0) then
	       nsumdtt=0
	       dt0=0.
	       mt0=0
	       do idp=1,nw
	         idr=nwi(idp)                      ! new
		     dt1=pdstime(id,idr)-pdst0(idr)    ! revised on 12/18/03
	         dt(idp)=dt1-dt0                   ! revised on 12/18/03
	         dt0=dt1                           ! revised on 12/18/03
               ndt(idp)=nint(dt(idp)/delta)
	         if(ndt(idp).lt.1) ndt(idp)=1
	         mt1=nint(pdstime(id,idr)/delta+0.5)    ! revised on 12/18/03
	         mt(idp)=mt1-mt0                        ! revised on 12/18/03
	         mt0=mt1                                ! revised on 12/18/03
	         move(idp)=int(mt(idp)/ndt(idp))
               nsumdtt=nsumdtt+ndt(idp)
	       enddo
	       nww=nw
             if(idr.lt.nre)then
	         nww=nww+1
		     dt(nww)=pdstime(id,nre)-pdst0(nre)-dt0
	         ndt(nww)=nint(dt(nww)/delta)
	         if(ndt(nww).lt.1) ndt(nww)=1           ! 04/20/04
	         mt(nww)=nint(pdstime(id,nre)/delta+0.5)-mt0
	         move(nww)=int(mt(nww)/ndt(nww))
               nsumdtt=nsumdtt+ndt(nww)
	       endif
	       lcc=0
	       mt0=0
	       do 120 k5=1,nww
               mm1=move(k5)     ! 1:mt1--every move(k5),remove delta (one point)
	         mt1=mt(k5)
	         mt2=mt0+mt1
	         do 110 k6=mt0+1,mt2
                 if(k6-mt0.ge.mm1) then
                   mm1=mm1+move(k5)
                   lcc=lcc+1
                 else
	             temp3(k6-lcc)=temp2(k6)
                 endif
110	         continue
	         mt0=mt2
120	       continue
             do 130 k8=mt0+1,leng1
               temp3(k8-lcc)=temp2(k8)
130	       continue
	       do 140 k9=leng1-lcc+1,leng1
	         temp3(k9)=0.0
140	       continue
	     else
	       nsumdtt=0
	       dt0=0.
	       mt0=0
	       do idp=1,nw
	         idr=nwi(idp)                      ! new
	         dt1=pdst0(idr)-pdstime(id,idr)    ! revised on 12/18/03
	         dt(idp)=dt1-dt0                   ! revised on 12/18/03
	         dt0=dt1                           ! revised on 12/18/03
               ndt(idp)=nint(dt(idp)/delta)
	         if(ndt(idp).lt.1) ndt(idp)=1
	         mt1=nint(pdstime(id,idr)/delta+0.5)    ! revised on 12/18/03
	         mt(idp)=mt1-mt0                        ! revised on 12/18/03
	         mt0=mt1                                ! revised on 12/18/03
	         move(idp)=int(mt(idp)/ndt(idp))
               nsumdtt=nsumdtt+ndt(idp)
	       enddo
	       nww=nw
             if(idr.lt.nre)then
	         nww=nww+1
		     dt(nww)=pdst0(nre)-pdstime(id,nre)-dt0
	         ndt(nww)=nint(dt(nww)/delta)
	         if(ndt(nww).lt.1) ndt(nww)=1           ! 04/20/04
	         mt(nww)=nint(pdstime(id,nre)/delta+0.5)-mt0
	         move(nww)=int(mt(nww)/ndt(nww))
               nsumdtt=nsumdtt+ndt(nww)
	       endif
	       lcc=0
	       mt0=0
	       do 160 k5=1,nww
               mm1=move(k5)     ! 1:mt1--every move(k5), add delta (one point)
	         mt1=mt(k5)
	         mt2=mt0+mt1
	         do 150 k6=mt0+1,mt2
                 if(k6-mt0.ge.mm1) then
                   mm1=mm1+move(k5)
                   lcc=lcc+1
		         temp3(k6+lcc-1)=temp2(k6)
  	             temp3(k6+lcc)=(temp2(k6-1)+temp2(k6)+temp2(k6+1))/3.
                 else
	             temp3(k6+lcc)=temp2(k6)
                 endif
150	         continue
	         mt0=mt2
160	       continue
             do 170 k8=mt0+1,leng1-lcc
               temp3(k8+lcc)=temp2(k8)
170	       continue
	     endif
         endif
180	 continue

	 if(idist.le.1)then
         do i5=1,leng1
           temp4(i5)=temp3(i5)
	   enddo
	   ltdt=0
	   leng4=leng1
	 else

C------------------------
c	 Do phase shift in time domain
         leng2=0
	   leng3=0                      ! new
	   ltdt=0
	   ldpoint=nint(depth(nre))
         lpdstime=nint(pdst0(nre)/delta)
	   lmargin=nint(dpdst0(nre)/delta)
	 do 190 i5=2,nre
         ltime=nint(pdst0(i5)/delta-pdst0(i5-1)/delta)
         ltd=nint(dpdst0(i5)/delta-dpdst0(i5-1)/delta)

        lsumt=0
	  leng4=leng3+ltime
	  if(ltd.ge.0) then
	  if(ltd.eq.0) then
	  do 205 i6=leng3+1,leng4
	  temp4(i6)=temp3(i6-leng3+leng2)
205	  continue
	  else	  
	  move3=nint(real(ltime)/real(ltd))
	  mm3=move3
	  do 200 i6=leng3+1,leng4
	  if((i6-leng3).ge.mm3) then
	  mm3=mm3+move3
	  lsumt=lsumt+1
	  else
	  temp4(i6-lsumt)=temp3(i6-leng3+leng2)
	  endif
200	continue
	  endif
	  leng2=leng2+ltime
	  leng3=leng3+ltime-ltd
	  ltdt=ltdt+ltd
         else
	   write(*,*) 'the program is error,stop here'
	   goto 220
         endif
190	continue
	leng4=leng1-ltdt
	do 230 i7=leng3+1,leng4
	temp4(i7)=temp3(i7-leng3+leng2)
230	continue

	endif             ! for if(idist.le.1)

	else
* idist < 0
         do i5=nbeg+1,npoint
           temp4(i5-nbeg)=temp1(i5)
	   enddo
	   ltdt=0
	   leng4=npoint-nbeg
	endif             ! for if(idist.ge.0)

* normalization
	if(inorm.ne.0)then
	  amp0=0.
	  do ii=1,leng4
	    amp=abs(temp4(ii))
	    if(amp.gt.amp0)amp0=amp
	  enddo
	  if(amp0.gt.1.e-12)amp0=1./amp0
	  do ii=1,leng4
	    temp4(ii)=temp4(ii)*amp0
	  enddo
	endif

	if(leng4.lt.np0)then
	  do ii=leng4+1,np0
	    temp4(ii)=0.
	  enddo
	endif
*      write(imgout,rec=ip) (temp4(ii),ii=1,np0)
      write(imgout,rec=ip) (temp4(ii)*am_cor,ii=1,np0)    ! 09/19/04
117   format(5f15.7)
	 close(12)
220	continue	 
	if(ip.eq.1)write(*,*)'ltdt =',ltdt,' leng1 =',leng1,
     &' leng4 =',leng4
	 return
	 end		
	 
