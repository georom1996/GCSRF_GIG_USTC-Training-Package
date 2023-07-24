*$ debug
c     ray tracing and calculating converted S piercing points
c     at chosen depth
C* Revised from pierc_new.f on 08/04/04
C* Using the ray parameter stored in each RF, not converted from gcarc
CCC************ZISTR**********
CCC	COMPUTING LOCATION OF CONVERSION POINTS FOR DIFFERENT DEPTH
cCC	FAR-FIELD FIRST-POINT APPROXIMATION
CCC	INPUT FROM SOURCE PARAMETER FILE
CCC	OUTPUT 
CCC	NW: THE NUMBER OF COMPUTING DEEPTH
CCC	J0: THE NUMBER OF LAYERS
CCC	DISTANCE=SUM(THICK*TG(VELOCITY*P0))
CCC*******************************
	include "stack.inc"
	parameter(nevtmax=1500,nstamax=500)            ! 01/02/05
      COMMON/RAYS/NA(nwmax,ICC),NRAY(nwmax,ICC),NEND,LMAX,NSP
      DIMENSION XO(nwmax),DEP(nwmax),TTO(nwmax),xox(nwmax),xoy(nwmax)
	dimension gcarc(num),stla(num),stlo(num),pp(num),baz(num)
	dimension evla(nevtmax),evlo(nevtmax)       !,evdp(num),az(num)
	dimension p0(icount),rang(icount)
      COMMON/TRAVEL/LP(ICC,2),LS(ICC,2) !07/26/10
      COMMON/STUFF/C(ICC),S(ICC),D(ICC),TH(ICC)
      COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC),NW,NWI(nwmax),
     &NWID(nwmax)
      CHARACTER*150 IAJ,IBJ,inname,inn1,inn2
      character*150 inname1(num),inname2(num)
	character*150 DIR,recfile,DIRSUB
      character*10 stfile
	integer ivar,igca,nst,ist,ldatadir,n0,nrf,ldirsub,nums,ldatadir0
	real var,varmin,varmax,gcarc0,amaxgca,sumgca,avergca,sumaz,averaz
	character*150 modfile,recordfile
	integer lrec,ndw(5)
	integer lensta(nstamax)               ! 01/02/05
	character*5 staname(nstamax)          ! 03/24/05
	integer numrf(nstamax)                ! 04/14/06
	
	! issue posted by Junliu Suwen >>>
	do i=1,nstamax
		lensta(i) = 0
	enddo
	! issue posted by Junliu Suwen <<<

	cont=pi/180.
	OPEN(5,FILE='pierc_new_n.in',STATUS='OLD')
*<1>
*	write(*,*)'output file name IAJ = ?'
	read(5,'(a)') ibj
	READ(5,500) IAJ
	WRITE(*,*) 'FILE:  ',IAJ
*<2> evla0, evlo0 
	read(5,'(a)') ibj
	read(5,*) evla0,evlo0
	write(*,*) 'the coordinate center of line',evla0,evlo0
*<3>
	read(5,'(a)') ibj
	read(5,*) np0,irayp
	irayp=0 
	write(*,*) 'output time point number =',np0
	write(*,*) 'ray parameter index =',irayp
C* irayp = 0: using user0 (pp(i)); else: using user1 (stored right after user0)
*<4>
*	WRITE(*,*) 'INPUT MODEL FILE NAME'
	read(5,'(a)') ibj
        READ(5,'(a)') modfile
        write(*,*)'model file =',modfile
500	FORMAT(A100)
510	FORMAT(1X,'FILE=',A12)
503	FORMAT(I2,4I3,F5.1,2F7.2,I5,F4.1,E10.4,3F7.2,F6.1)
*<5>
	read(5,'(a)')ibj
	read(5,*)ivar,varmin,varmax
	write(*,*)'ivar=',ivar,'  range:',varmin,' --',varmax
*<6>
C	WRITE(*,*) 'READ NW,(NWI(I),NWID(I),I=1,NW)
	read(5,'(a)') ibj
* ray number, ray indexes, corresponding depth indexes
	READ(5,*) NW,(NWI(I),NWID(I),I=1,NW)   
*<7>
C	WRITE(*,*) 'READ five depth nw indexs for outputting piercing points'
	read(5,'(a)') ibj
* indexs in nwi: ndw(1:5)
	READ(5,*) (ndw(i),i=1,5)   

	nw0=nw/2+1
	rlat=1./rad
*	rlon=1./(rad*cos(evla0*cont))             ! changed on 10/28/04
	
C***** modfile: MODEL FILE*******
      OPEN(9,FILE=modfile,status='old')
      READ(9,*)JO
      DO J=1,JO
      READ (9,*) C(J),S(J),D(J),TH(J)
      END DO
      CLOSE(9)
      DO 851 J=1,JO
      CC(J)=C(J)
      SS(J)=S(J)
      DD(J)=D(J)
851   TTH(J)=TH(J)

	YY=0.
	DO K=1,NW
	IN=NWID(K)
	I1=1
	IF(K.GT.1) I1=NWID(K-1)+1
	DO I=I1,IN
	YY=YY+TH(I)
	END DO
	DEP(K)=YY
C     READ(5,*) (NA(K,J),J=1,JO)
	END DO


      LMAX = JO
      CALL CURAY(JO)
      MODP=5
	MODS=3
	NENN=JO
      NEND=NENN
      LFINAL=NW

C***** '' *******
*<8>
C	WRITE(*,*) 'READ directory containing RFs'
	read(5,'(a)') ibj
	read(5,'(a)')DIR
	ldir=lengths(DIR)
*<9>
C	WRITE(*,*) 'READ number of subdirectories'
	read(5,'(a)') ibj
	read(5,*) nst0                                     ! 04/21/06
	nst=iabs(nst0)                                     ! 04/21/06
	write(*,*)'subdirectory number =',nst

        sumgca=0.
	sump=0.
	aminp=10.
        sumaz=0.
	n0=0
	n=0
	nums=0
	ldatadir0=0
	nevt=0
*04/14/06
	nums0=0
	do i=1,nstamax
	  numrf(i)=0
	enddo
* end
	do 1000 ist=1,nst
	read(5,'(a)')DIRSUB
	read(5,'(a)')recordfile
	ldirsub=lengths(DIRSUB)
	lrec=lengths(recordfile)
	
	recfile=DIR(1:ldir)//DIRSUB(1:ldirsub)//'/'//recordfile(1:lrec)
	write(*,'(a)')recfile
	ldir2=ldir+ldirsub+1
	open(10,file=recfile,status='old')
	do 300 istt=1,10000
	read(10,'(a)')stfile
	write(*,'(a)')stfile
      if(stfile(1:2).eq.'  ') goto 99
	ldatadir=lengths(stfile)
	nums0=nums+istt                                ! 04/14/06
	lensta(nums0)=ldir2                            ! 01/02/05
	staname(nums+istt)=stfile(1:ldatadir)          ! 03/24/05
	if(ldatadir0.lt.ldatadir)ldatadir0=ldatadir
	read(10,*)nrf
	n0=n0+nrf
      do 310 i=1,nrf
      read(10,111,end=99) inname
111   format(a100)
	linname=lengths(inname)
      if(inname(1:2).eq.'  ') then
	  n0=n0-1
	  goto 310
	endif
	n1=n+1
	nevt1=nevt+1
      inname1(n1)=recfile(1:ldir2)//stfile(1:ldatadir)//'/'//
     &inname(1:linname-4)//'.bhz'
      inname2(n1)=recfile(1:ldir2)//stfile(1:ldatadir)//'/'//
     &inname(1:linname)
	lengthr=5
	nrecl=4*lengthr
C* 04/21/06
	if(nst0.lt.0) then
      open(12,file=inname2(n1),status='old',recl=nrecl,err=9,
     &	access='direct',convert='BIG_ENDIAN')
	else
      open(12,file=inname2(n1),status='old',recl=nrecl,err=9,
     &	access='direct')
	endif
C* end 04/21/06
	read(12,rec=7) test1,stla(n1),stlo(n1),stel,stdp
	read(12,rec=8) evla(nevt1),evlo(nevt1),evel,evdp,test2
	read(12,rec=9) pp(n1),temp,test1,test2,sss1
	if(irayp.ne.0) then
	pp(n1)=temp
	else
	pp(n1)=pp(n1)
	endif
	
      read(12,rec=11) dis1,az,baz(n1),gcarc(n1),sss1
      close(12)
	if(ivar.eq.0)then
	  var=dis1
	elseif(ivar.eq.1)then
	  var=gcarc(n1)
	else
	  var=baz(n1)
	  if(var.lt.varmin.and.varmax.gt.360.)then
          var=var+360.
	  endif
	endif
	if(var.ge.varmin.and.var.le.varmax)then
          n=n1
	  numrf(nums0)=numrf(nums0)+1           ! 04/14/06
C          sumgca=sumgca+gcarc(n)
          if(pp(n).lt.aminp)aminp=pp(n)
	  sump=sump+pp(n)
          sumaz=sumaz+baz(n)
          sumgca=sumgca+gcarc(n)
113   format(5G15.7)
112   format(a80)
        
        

C* 08/20/04
	do j=1,nevt
	if(abs(evla(nevt1)-evla(j)).lt.0.01.and.
     & abs(evlo(nevt1)-evlo(j)).lt.0.01) goto 310
	enddo
	nevt=nevt1
C* end 08/20/04

	endif
	goto 310
9	write(*,'(a)')inname2(n1)
310   continue
	read(10,'(a)')
300	continue
99    close(10)
	nums=nums+istt-1
1000	continue
	close(5)
        avergca=sumgca/n
	averp=sump/n
        averaz=sumaz/n
	ldatadir0=min(ldatadir0,5)

	pp0=averp
	gcarc0=avergca
      
	write(*,*)'Total RF number =',n0,'  used RF number =',n

        write(*,*) '***************'
	OPEN(7,FILE=IAJ,STATUS='UNKNOWN')
	open(17,file='sta_name.dat',status='unknown')
	open(18,file='station.dat',status='unknown')
	open(19,file='depth1.dat',status='unknown')
	open(20,file='depth2.dat',status='unknown')
	open(21,file='depth3.dat',status='unknown')
	open(22,file='depth4.dat',status='unknown')
	open(23,file='depth5.dat',status='unknown')
C* 08/20/04
	open(24,file='events.txt',status='unknown')
	do i=1,nevt
	  write(24,115)evlo(i),evla(i)
	enddo
115	format(f10.4,2x,f10.4)
	close(24)
C* 08/20/04

	WRITE(7,*) NW,(NWI(I),DEP(I),I=1,NW)
	write(7,*) evla0,evlo0
	write(7,*) gcarc0,pp0,aminp,averaz	

	write(7,*) n, np0
	write(7,'(a)')DIR(1:ldir)
      inn1=' '
	JJ=0
	nums1=1                              ! 01/02/05
	do 897 ii=1,n
	inn2=inname1(ii)
	dx=0.
	dy=0.
C* get distance (in km) between receiver and reference point in north and east direc.
	call guass(stla(ii),stlo(ii),evla0,evlo0,dx,dy)
	rlon=1./(rad*cos(stla(ii)*cont))                   ! 10/28/04
	
	inname=inname2(ii)
	linname=lengths(inname)
	if(linname.gt.149)then
	  write(*,*)'The file name is too long: ',inname(ldir+1:linname)
	  goto 1001
	endif
      WRITE(7,'(a)') inname(ldir+1:linname)
      WRITE(7,*) ii,dx,dy
	ldir2=lensta(nums1)                  ! 01/02/05
      if(inn1(ldir2+1:ldir2+5).ne.inn2(ldir2+1:ldir2+5)) then
	  inn1=inn2
	  nums1=nums1+1                      ! 01/02/05
	  write(*,*)'nums1 =',nums1          ! 01/02/05
	  write(18,*) stla(ii),stlo(ii),dx,dy
	  write(17,124) staname(nums1-1)
	endif
114	format(1x,2f10.4,3x,a5)
124	format(3x,a5)
        ppo=pp(ii)
        range=gcarc(ii)
        azb=baz(ii)
        JJ=JJ+1
2	CONTINUE
	NNR=JJ 
      DO NSP=1,LFINAL	   
	DO J=1,NENN
CL 07/26/10
c      LP(J)=0
c	LS(J)=0
c	IF(J.gt.NWID(NSP))then
c	LP(J)=1
c	ELSE
c	LS(J)=1
c	END IF
        LP(J,1)=0
        LP(J,2)=0
        LS(J,1)=0
        LS(J,2)=0
        IF(J.le.NWID(NSP))then
            LP(J,1)=-1
            LS(J,2)=1
	ENDIF
CL 07/26/10 end
        END DO
      X=RANGE*rad
      PO=PPO
      CALL FIND (JO,PO,XO)        ! get offset from the convert point at surface
      CALL FIND2 (JO,PO,TTO)      ! travel time
      END DO
C* xox: offset of convert point to reference point in north
C* xoy: offset of convert point to reference point in east
	do 893 i2=1,nw
	xox(i2)=xo(i2)*cos(azb*cont)+dx
	xoy(i2)=xo(i2)*sin(azb*cont)+dy
893	continue
	acl1=cos(azb*cont)
	acl2=sin(azb*cont)
	do lc=1,5
	  write(18+lc,*) xo(ndw(lc))*acl1*rlat+stla(ii),
     & xo(ndw(lc))*acl2*rlon+stlo(ii),xox(ndw(lc)),xoy(ndw(lc))
	enddo

      WRITE(7,512) RANGE,azb,cos(azb*cont),PO
	WRITE(7,511) (xox(I),I=1,NW),(xoy(i),i=1,nw)
	write(7,511) (TTO(I),I=1,NW)
897   continue
1002	CONTINUE
1001	CONTINUE
	WRITE(*,*) 'NNR',NNR
	write(*,*)'total number of stations =',nums,nums1      ! 01/02/05
	write(*,*)'output piercing points at 5 depths:'
	write(*,*)(DEP(ndw(I)),I=1,5)
	write(*,*)'total event number =',nevt
511	FORMAT(1X,10F8.2)
512	FORMAT(1X,2F8.2,2E12.4)
      CLOSE(5)
	close(21)
	close(20)
	close(19)
	close(18)
	close(17)
	close(22)
	close(23)
	CLOSE(7)

      STOP
      END
