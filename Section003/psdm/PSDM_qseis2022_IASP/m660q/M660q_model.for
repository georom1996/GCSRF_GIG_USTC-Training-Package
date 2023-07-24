*##########################################################*
*                                                          *
*  M660q ---- calculate the travel time of free-surface    *
*             reflections or discontinuity conversions     *
*             for each specified depth                     *
*                                                          *
*##########################################################*
C* Using cwbq model as the reference model
C* Ray parameters for different epi. distances (nmo.dat) are obtained
C* based on either PREM or iasp91 model. It seems no much difference from IASP91.

*$debug
       parameter(icc=500,ina1=150,lgreen=1050)
	 parameter(depmax=1200,rad=111.195)
       common/rays/na(ina1),nray(ina1),nend,lmax,nsp,love,
     *   nzr,iwh,jo
       common/travel/lp(icc),ls(icc),ll,fdp,map,nplnw
       common/stuff/c(icc),s(icc),d(icc),th(icc)
       common/orst/cc(icc),ss(icc),dd(icc),tth(icc)
       common/mod/cka2(icc),ckb2(icc),cmu(icc),roh(icc)	 
       common/save/rd(2,2,icc),td(2,2,icc),ru(2,2,icc),tu(2,2,icc)
     *  ,a(icc),b(icc),term(2,2),rsur(2,2)
	 dimension range(150),ppo(150),x(icc),depth(icc),dep(icc)
	 dimension na0(ina1,icc),nray0(ina1,icc),ilay0(ina1)
	 dimension tuamp(ina1,icc),time(ina1,icc),syn(ina1,lgreen)
	 dimension dtime(ina1,icc),dv(icc),dv3(icc)
       complex rd,td,ru,tu,ckb2,cka2
       complex term,rsur,a,b,tuamp,tutemp,tulay
       character*100 velmod,outfile,rayfile,chartemp
       integer iflat,itype

	open(1,file='m660q_model.in',status='old')
	read(1,'(a)')chartemp
	read(1,'(a)')velmod
	write(*,*)'input velocity model file: ',velmod
	read(1,'(a)')chartemp
	read(1,'(a)')rayfile
	write(*,*)'input ray file: ',rayfile
	read(1,'(a)')chartemp
	read(1,'(a)')outfile
	write(*,*)'output traveltime file: ',outfile
	read(1,'(a)')chartemp
	read(1,*) iflat,itype
CL 07/27/10: itype = 0: free-surface reflection; else: conversion (>0: Ps; <0: Sp)
CL 05/04/12: cannot be applied to the itype < 0 case
	write(*,*)'flattening index iflat (=0: no flattening; else: do): ',iflat
        write(*,*)'phase type: ',itype
	close(1)

       open(12,file=outfile,status='unknown')
       open(9,file=velmod,status='old')
       read(9,*)joo
c	write(12,*) joo
	 depth0=0.
       do 11 j=1,joo
        read (9,*) c(j),s(j),d(j),th(j)      ! vp,vs,density,depth interval
	  depth(j)=depth0+th(j)
	  if(depth0.ge.depmax) goto 13         ! depmax instead of "660"
	  depth0=depth0+th(j)
11	 continue
13	 jo=j                                  ! jo is the index for the maxi. dep
	 close(9)
	
       open(9,file='nmo.dat')
       read (9,*) nrange
       do i=1,nrange
       read (9,*)range(i),ppo(i)       ! epi. distance, ray parameter
       enddo
       close(9)
* new
	 nrange=nrange+1
	 range(nrange)=180.
	 ppo(nrange)=0.
	 
      do 21 j=1,jo
       cc(j)=c(j)
       ss(j)=s(j)
       dd(j)=d(j)
       tth(j)=th(j)
21	continue	
      lmax = jo
      if(iflat.ne.0) call curay(jo)
* Ling 08/11/03
	dv0=0.
	dv1=0.
	do j=1,jo
	  dv(j)=dv0+(cc(j)-ss(j))*th(j)
	  dv0=dv(j)
* Ling 08/22/03
	  dv3(j)=dv1+(cc(j)**3-ss(j)**3)*th(j)
	  dv1=dv3(j)
*wrong	  dv0=dv0+(1./ss(j)-1./cc(j))*th(j)
	enddo
* end
	
*	open(10,file='mray_mod.dat')
	open(10,file=rayfile,status='old')
	read(10,*)nre
*	write(12,*) nre,nrange
	do 31 iray=1,nre
	 read(10,32)ilay0(iray),(na0(j,iray),j=1,ilay0(iray))
	 read(10,33)(nray0(j,iray),j=1,ilay0(iray))
C* 07/12/10: for free-surface reflected phases:
         if(itype.eq.0)then                        ! CL 07/27/10
	   dep(iray)=0.
           do j=ilay0(iray),2,-1
             dep(iray)=dep(iray)+th(na0(j,iray))
             if(na0(j,iray).eq.na0(j-1,iray)) goto 31
           enddo
           if(iray.eq.1) dep(iray)=0.              ! CL 06/15/11
CL: I find before this date, dep(1) calculated is the maximum depth, not 0.
CL: This doesn't affect my results (may cause errors if call the subroutine moveout_d), but should be corrected
         endif
31	continue
32	format(i2,1x,100i2)
33	format(3x,100i2)	
	close(10)
	
C*	open(13,file='m660q.dat',status='unknown')	
	open(14,file='m660q.rst',status='unknown')
C*	write(13,*)nrange
	dt=0.1
	am=10
	nl=2**am	
	tstart=5
	nsta=tstart/dt
		
	modp=5
	mods=3

	do 41 ir=1,nrange	
	do 51 iray=1,nre
	 jotemp=ilay0(iray)
	 do 53 ijo=1,jo
	  lp(ijo)=0
	  ls(ijo)=0
	  lcount=0
53	 continue	 
	 do 61 j=1,jotemp
	  na(j)=na0(j,iray)
	  jt=na(j)
	  nray(j)=nray0(j,iray)
	  if(nray(j).eq.modp)then
	   lp(jt)=lp(jt)+1	     
	  else
	   ls(jt)=ls(jt)+1
	  end if 
	  if(itype.gt.0)then                       ! CL 07/27/10
            lcount=lcount+ls(jt)
          else
            lcount=lcount+lp(jt)
          endif                                    ! CL 07/27/10
61	 continue
c          if(ir.eq.1) write(*,*) lcount
C* for Ps phases
         if(itype.ne.0)then                        ! CL 07/27/10
	   dep(iray)=0.
	   if(lcount.ne.0) then
	     do 63 k5=1,lcount
               dep(iray)=dep(iray)+th(k5)     ! depth ?
63	     continue
	   endif
         endif                                     ! CL 07/27/10
	 po=ppo(ir)/rad
       call find (po,to,x)
	 time(ir,iray)=to                ! to -- trave time
*	 write(*,*)'time =',to

      po2=po*po
      pa=3.1415926
* Ling 08/11/03
	 dtime(ir,iray)=0.
* Ling 08/22/03: second part 
	 if(lcount.ne.0) dtime(ir,iray)=
     & (dv(lcount)+dv3(lcount)*0.25*po2)*po2*0.5
*wrong	 if(lcount.ne.0)dtime(ir,iray)=to-dv(lcount)

      do 81 i=1,jo
       cka2(i)=1.0/cc(i)/cc(i)         ! 1./vp.^2
       ckb2(i)=1.0/ss(i)/ss(i)         ! 1./vs.^2
       cmu(i)=ss(i)*ss(i)*dd(i)        ! density*vs.^2
81    continue
      
*      call reftra(po,po2,jo)	         ! reflection & transmision coef.
*	tutemp=(1.,0.)
*	do 101 i=2,jotemp
*	 nlay=na(i)+1
*	 mod1=nray(i-1)
*	 mod2=nray(i)
*	 if(na(i).lt.na(i-1))then
*	  if(mod1.eq.modp.and.mod2.eq.modp)tulay=tu(1,1,nlay)
*	  if(mod1.eq.modp.and.mod2.eq.mods)tulay=tu(1,2,nlay)
*	  if(mod1.eq.mods.and.mod2.eq.modp)tulay=tu(2,1,nlay)
*	  if(mod1.eq.mods.and.mod2.eq.mods)tulay=tu(2,2,nlay)
*	 else if(na(i).gt.na(i-1))then
*	  if(mod1.eq.modp.and.mod2.eq.modp)tulay=td(1,1,nlay)
*	  if(mod1.eq.modp.and.mod2.eq.mods)tulay=td(1,2,nlay)
*	  if(mod1.eq.mods.and.mod2.eq.modp)tulay=td(2,1,nlay)
*	  if(mod1.eq.mods.and.mod2.eq.mods)tulay=td(2,2,nlay) 	 
*	 else if(na(i).eq.na(i-1).and.na(i).eq.1)then
*	  if(i/2*2.ne.i)then
*	   if(mod1.eq.modp.and.mod2.eq.modp)tulay=rsur(1,1)
*	   if(mod1.eq.modp.and.mod2.eq.mods)tulay=rsur(1,2)
*	   if(mod1.eq.mods.and.mod2.eq.modp)tulay=rsur(2,1)
*	   if(mod1.eq.mods.and.mod2.eq.mods)tulay=rsur(2,2) 	 
*	  else
*	   if(mod1.eq.modp.and.mod2.eq.modp)tulay=rd(1,1,nlay)
*	   if(mod1.eq.modp.and.mod2.eq.mods)tulay=rd(1,2,nlay)
*	   if(mod1.eq.mods.and.mod2.eq.modp)tulay=rd(2,1,nlay)
*	   if(mod1.eq.mods.and.mod2.eq.mods)tulay=rd(2,2,nlay) 	  
*	  end if	  
*	 else if(na(i).eq.na(i-1).and.na(i).ne.1)then
*	  if(mod1.eq.modp.and.mod2.eq.modp)tulay=rd(1,1,nlay)
*	  if(mod1.eq.modp.and.mod2.eq.mods)tulay=rd(1,2,nlay)
*	  if(mod1.eq.mods.and.mod2.eq.modp)tulay=rd(2,1,nlay)
*	  if(mod1.eq.mods.and.mod2.eq.mods)tulay=rd(2,2,nlay) 	  
*	 end if
*	 tutemp=tutemp*tulay
*101	continue
*c	for Rpr and Rsr	
*	if(nray(jotemp).eq.modp)tutemp=tutemp*term(1,1)
*	if(nray(jotemp).eq.mods)tutemp=tutemp*term(1,2)
*c	for Rpz and Rsz	
*c	if(nray(jotemp).eq.modp)tutemp=tutemp*term(2,1)
*c	if(nray(jotemp).eq.mods)tutemp=tutemp*term(2,2)	
*	tuamp(ir,iray)=tutemp               ! amplitude
51	continue
	if(ir.eq.1) then
	write(12,1001) nre,nrange
	write(12,1002) (dep(k6),k6=1,nre)
	else
	endif	
	write(12,1003) range(ir),po
	write(12,1004)(time(ir,i)-time(ir,1),i=1,nre)
C*	write(12,1004)(dtime(ir,i),i=1,nre)
1001	format(2i6)
1002	format(5f15.3)
1003	format(f15.2,f15.7)
1004	format(5f15.5)
*c	write(12,*)(tuamp(ir,i),i=1,nre)
*	do 105 i=1,nre
*	if(abs(real(tuamp(ir,i))).ge.1e-02)then
*         if(ilay0(i).gt.42) goto 1071
c         write(14,*) 'range(ir),po=',range(ir),po
         write(14,*) range(ir),po
         write(14,*) time(ir,i)
c         write(14,32)ilay0(i),(na0(j,i),j=1,ilay0(i))
c         write(14,33)(nray0(j,i),j=1,ilay0(i))
c         write(14,*)'time,amp=',time(ir,i),tuamp(ir,i)
*	endif
1071  continue        
105	continue
*	do 111 it=1,nl		
*	 syn(ir,it)=0.0
*111	continue
*	 syn(ir,nsta)=real(tuamp(ir,1))
*	 do 121 iray=2,nre	 
*	  tint=time(ir,iray)-time(ir,1)+tstart
*	  ntint=tint/dt
*	  if(ntint.le.nl)
*     & syn(ir,ntint)=syn(ir,ntint)+real(tuamp(ir,iray))	  
*121	 continue	
*	write(13,*)nl,dt
*	write(13,*)(syn(ir,i),i=1,nl) 	 
*	if(ir.eq.1)then
*	open(15,file='m660q0.dat')
*	write(15,*)nl
*	do inl=1,nl
*	write(15,*)inl*dt-tstart,syn(ir,inl)/0.1-0.08
*	end do
*	close(15)
*	end if			
41	continue	
	close(12)
C*	close(13)
        close(14)
 	
C*	write(*,*)'output file:m660q.dat'
	write(*,*)'output file:',outfile    
	write(*,*)'output file:m660q.rst'	
	stop
	end


      SUBROUTINE CURAY(JO)
      PARAMETER(ICC=500)
      COMMON/STUFF/C(ICC),S(ICC),D(ICC),TH(ICC)
      COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC)
      COMMON/DEPTHS/DEPTH(ICC)

* calculate the depth of middle point for each layer
      DEPTH(1)=TTH(1)/2.0
      DO 10 J = 2,JO
  	 DEPTH(J)=DEPTH(J-1)+(TTH(J)+TTH(J-1))/2.
10	CONTINUE
      DO 5 J = 1,JO
       Q = 6371.0 / (6371.0-DEPTH(J))
       CC(J)=(C(J))*Q              ! vp*Q ?
       SS(J)=(S(J))*Q              ! vs*Q ?
       DD(J)=(D(J))*Q              ! density*Q ?
       TTH(J)=(TH(J))*Q            ! depth interval * Q ?
5     CONTINUE
     
      RETURN
      END

      SUBROUTINE FIND (PO,TO,X)
      PARAMETER(ICC=500,INA1=150)
      COMMON/RAYS/NA(INA1),NRAY(INA1),NEND,LMAX,NSP,LOVE,
     *   NZR,IWH,JO
      COMMON/TRAVEL/LP(ICC),LS(ICC),LL,FDP,MAP,NPLNW
      COMMON/STUFF/C(ICC),S(ICC),D(ICC),TH(ICC)
      COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC)
      COMMON/TSTAR/TQ,TVL,QA(ICC),QB(ICC)
      REAL*8 E(200),F(200),XX,YY,P
	DIMENSION X(ICC)
      P = PO
      DO 10 J = 1,JO
	 XX=DABS(1./SS(J)-P)
	 YY=DABS(1./SS(J)+P)
       F(J) =DSQRT(XX*YY)             ! sqrt(1/vs^2-p^2)
	 YY=DABS(1./CC(J)+P)
	 XX=DABS(1./CC(J)-P)
       E(J) =DSQRT(XX*YY)             ! sqrt(1/vp^2-p^2)
 10   CONTINUE
      TOTEM = 0.0
	X0=0.0
	DO 11 J=1,JO-1
	X0=X0+(TTH(J)*LP(J)/E(J)+TTH(J)*LS(J)/F(J))*PO
	X(J)=X0
	TOTEM=TOTEM+E(J)*TTH(J)*LP(J)+F(J)*TTH(J)*LS(J)
11	CONTINUE
	TO=TOTEM	
      RETURN
      END

      SUBROUTINE REFTRA(XK,XK2,N)
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=500)
      COMPLEX RD,TD,RU,TU,CKB2,CKA2,CR2,A1,A2,B1,B2,AB1,AB2,AB12,Z1
      COMPLEX ROH1,ROH2,CD
      COMPLEX A1B2,B1A2
      COMPLEX TERM,RSUR,Z,CQ1,CQ2,CQ3,CQ4,Q1,Q2,Q3,Q4,Q5,Q6,D,OMEGA
      COMPLEX ROH12,A,B
      COMMON/ORST/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL)
      COMMON/MOD/CKA2(NLL),CKB2(NLL),CMU(NLL),ROH(NLL)
      COMMON/SAVE/RD(2,2,NLL),TD(2,2,NLL),RU(2,2,NLL),TU(2,2,NLL)
     $,A(NLL),B(NLL),TERM(2,2),RSUR(2,2)
C
      Z=CKA2(1)
      A(1)=CR2(XK2,Z)       ! i*cos(ang_P1)
      Z=CKB2(1)
      B(1)=CR2(XK2,Z)       ! i*cos(ang_S1)
      DO 100 J=1,N
      IF(J.EQ.1) GOTO 199
      Z=CKA2(J)
      A(J)=CR2(XK2,Z)       ! i*cos(ang_Pj)
      Z=CKB2(J)
      B(J)=CR2(XK2,Z)       ! i*cos(ang_Sj)
      I1=J-1
      I2=J
      A1=A(I1)
      B1=B(I1)
      A2=A(I2)
      B2=B(I2)
      IF(ABS(AFA(I2)-AFA(I1)).LT.0.001) GOTO 200
      AB1=A1*B1
      ROH2=CKB2(I2)*CMU(I2)    ! density
      AB2=A2*B2
      AB12=AB1*AB2
      CQ1=CMU(I2)-CMU(I1)
      ROH1=CKB2(I1)*CMU(I1)
      ROH12=0.5*ROH1
      CQ4=-CQ1*XK2-ROH12
      Z=0.5*ROH2
      CQ2=CQ4+Z
      CQ3=CQ1*XK2-Z
      Q1=XK2*AB12*CQ1*CQ1
      Q2=XK2*CQ2*CQ2
      Q3=AB1*CQ3*CQ3
      Q4=AB2*CQ4*CQ4
      ROH12=ROH12*Z
      A1B2=A1*B2
      Q5=ROH12*A1B2
      B1A2=B1*A2
      Q6=ROH12*B1A2
      D=Q1+Q2-Q3-Q4-Q5-Q6
C
      X=CABS(D)
      IF(X.GT.1.0E-15) GOTO 606
      ISAFE=1
      WRITE(*,605)
605   FORMAT(1X,'D=0.0! STOP IN REFTRA')
      STOP
606   CONTINUE
C
      Z=CQ3*B1+CQ4*B2
      Z=CD(Z,D)
      TD(1,1,I2)=ROH1*A1*Z
      TU(1,1,I2)=ROH2*A2*Z
      Z=XK*(CQ1*A1B2+CQ2)
      Z=CD(Z,D)
      TD(1,2,I2)=ROH1*B1*Z
      TU(2,1,I2)=ROH2*A2*Z
      Z=XK*(CQ1*B1A2+CQ2)
      Z=CD(Z,D)
      TD(2,1,I2)=ROH1*A1*Z
      TU(1,2,I2)=ROH2*B2*Z
      Z=A1*CQ3+CQ4*A2
      Z=CD(Z,D)
      TD(2,2,I2)=ROH1*B1*Z
      TU(2,2,I2)=ROH2*B2*Z
      Z=Q1-Q2-Q3+Q4
      Z1=Z-Q5+Q6
      RD(1,1,I2)=CD(Z1,D)
      Z1=Z+Q5-Q6
      RD(2,2,I2)=CD(Z1,D)
      Z=Q1-Q2+Q3-Q4
      Z1=Z+Q5-Q6
      RU(1,1,I2)=CD(Z1,D)
      Z1=Z-Q5+Q6
      RU(2,2,I2)=CD(Z1,D)
      Z=-2.*XK*(CQ3*CQ2-AB2*CQ1*CQ4)
      Z=CD(Z,D)
      RD(1,2,I2)=B1*Z
      RD(2,1,I2)=A1*Z
      Z=-2.*XK*(CQ4*CQ2-AB1*CQ1*CQ3)
      Z=CD(Z,D)
      RU(1,2,I2)=B2*Z
      RU(2,1,I2)=A2*Z
      GOTO 100
C
199   CONTINUE
      A1=A(1)
      B1=B(1)
      AB1=A1*B1
      OMEGA=XK2-0.5*CKB2(1)
      Z=XK2*AB1-OMEGA*OMEGA
C
      X=CABS(Z)
      IF(X.GT.1.0E-20) GOTO 906
      ISAFE=1
      WRITE(*,907)
907   FORMAT(1X,'Z=0.0! STOP IN REFTRA FOR FIRST LAYER')
      STOP
906   CONTINUE
      Z1=CKB2(1)*XK*AB1
      Z1=CD(Z1,Z)
      TERM(1,1)=Z1
      TERM(2,2)=Z1
      Z1=CKB2(1)*OMEGA
      Z1=CD(Z1,Z)
      TERM(1,2)=B1*Z1
      TERM(2,1)=A1*Z1
      Z1=XK2*AB1+OMEGA*OMEGA
      Z1=CD(Z1,Z)
      RSUR(1,1)=Z1
      RSUR(2,2)=Z1
      Z1=2.*XK*OMEGA
      Z1=CD(Z1,Z)
      RSUR(1,2)=B1*Z1
      RSUR(2,1)=A1*Z1
      GOTO 100
200    Z=CMPLX(0.,0.)
      Z1=CMPLX(1.,0.)
      DO 80 K1=1,2
      DO 80 K2=1,2
      RD(K1,K2,I2)=Z
      RU(K1,K2,I2)=Z
      TD(K1,K2,I2)=Z
80    TU(K1,K2,I2)=Z
      TD(1,1,I2)=Z1
      TD(2,2,I2)=Z1
      TU(1,1,I2)=Z1
      TU(2,2,I2)=Z1
100   CONTINUE
      RETURN
      END

      COMPLEX FUNCTION CR2(P2,C2)
      DOUBLE PRECISION U,X,R,W1,W2,R1,R2
      COMPLEX CZ,C2
      CZ=C2-P2
      U=REAL(CZ)
      X=AIMAG(CZ)
      R=DSQRT(X*X+U*U)
      W1=DABS(R+U)/2.
      W2=DABS(R-U)/2.
      R1=DSQRT(W1)
      R2=DSQRT(W2)
      R11=R1
      R22=R2
      CR2=R22+R11*CMPLX(0.,1.)
      RETURN
      END

       COMPLEX FUNCTION CD(Z1,Z2)
       COMPLEX Z1,Z2
       DOUBLE PRECISION X1,X2,Y1,Y2,W1,W2,W3
       X1=REAL(Z1)
       Y1=AIMAG(Z1)
       X2=REAL(Z2)
       Y2=AIMAG(Z2)
       W3=X2*X2+Y2*Y2
       W1=(X1*X2+Y1*Y2)/W3
       W2=(Y1*X2-Y2*X1)/W3
       W11=W1
       W22=W2
       CD=W11+W22*(0.,1.)
       RETURN
       END
