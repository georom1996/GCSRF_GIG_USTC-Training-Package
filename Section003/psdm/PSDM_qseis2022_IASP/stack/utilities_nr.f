       subroutine findrec(x1,y1,inumb,ninw,inw0)
	 include "stack.inc"
	 common/param1/baz(num),gcarc(num),stla(num),stlo(num)
       common/param2/ppo(num),icnt(num),n,ldir
	 common/param3/piercx(num,nwmax),piercy(num,nwmax),tto(num)
       COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC),NW,NWI(nwmax),
     & NWID(nwmax)
	 common/locate1/begla,beglo,endla,endlo,stepbin,SA,CA,leasttra,
     & XBIN(nwmax),YBIN0(nwmax),DYBIN(nwmax),YBINM(nwmax)
	 common/locate2/numsta,numtra,numbin,dep(nwmax),
     & avergca,averaz	 
       common/param4/nsta(ina,nwmax),ntra(ina,nwmax),depth(icc),
     & YBIN(ina,nwmax)
	 character*150 iaj,ibj,inn1,inn2
	 integer lbiname(maxnum)
	 integer ninw,inw0(nwmax)

	 cont=pi/180.
	 j1=(inumb-1)*ninw
       do 1000 jcl=1,ninw
	   j=inw0(jcl)
         icon=0
	   ybin1=0.
	   ybin2=YBIN0(jcl)
	   do 1005 iy=1,1000
         do 1010 k=1,n
	   px=piercx(k,j)
	   py=piercy(k,j)
         disX=abs((pX-X1)*ca-(pY-Y1)*sa)
         disY=abs((pY-Y1)*ca+(pX-X1)*sa)
         if(disX.le.XBIN(jcl).AND.DISY.Lt.ybin2.and.DISY.GE.ybin1) then
            icon=icon+1
	     lbiname(icon)=k                     ! new
         endif
1010    continue
	  if(icon.ge.numtra.or.ybin2.ge.YBINM(jcl)) goto 1020
	  ybin1=ybin2
	  ybin2=ybin2+DYBIN(jcl)
1005	  continue
1020    write(8,*) 'trace number=',icon,' YBIN =',ybin2
        ntra(inumb,jcl)=icon
	  YBIN(inumb,jcl)=ybin2
	  j1=j1+1
	  write(10,rec=j1)(lbiname(j2),j2=1,maxnum)   ! new
1000   continue
	                              
	 return
         end


c       ***************************
c       get string's actually length
c       ***************************
        Integer Function Lengths(str)
        Character*(*) str
        do 15 i = len(str),1,-1
                if(str(i:i) .NE. ' ') go to 20
15      continue
20      Lengths = i                 ! modified
        End


c----------------------------------------
c     convert an integer to a string
c----------------------------------------
      subroutine num2str(inum,str,len)
c-------------------------------------------------------------------
c Input
c     inum    --- integer num
c Output
c     str     --- converted string
c     len     --- effective length of the string
c-------------------------------------------------------------------
c     Written by Ling Chen at UCSC on September, 2001
c

      parameter (nemax=10)
      integer inum,len
      character str(nemax+1)
      integer i,j,k,m,n,inum1
      
      inum1=inum
      k=0
      if(inum1.lt.0)then
        str(1)='-'
        k=k+1
        inum1=inum1*(-1)
      endif
      do i=1,nemax
        n=10**i
        if(inum1.lt.n)goto 10
      enddo

10    len=k+i
      do j=i-1,0,-1
        n=10**j
        m=inum1/n
        k=k+1
        str(k)=char(m+48)
        inum1=inum1-m*n
      enddo

      return
      end


	SUBROUTINE GUASS(FA,LA,FA0,LA0,X,Y)
	parameter(pi=3.1415926536,nll=15)
      REAL*4 A,B,CONT,E,E2,EP,EP2
	REAL LA,LA0,L,N,L2,L3,L4,L5,L6
	a=6378.1370000
	b=6356.7523142
      cont=pi/180.0
      SM=6367.559*PI/2.0
      SM=SM*(FA-FA0)/90.0
	E2=(A*A-B*B)/A/A
	EP2=(A*A-B*B)/B/B
	EP=SQRT(EP2)
	E=SQRT(E2)
	L=(LA-LA0)*CONT
	X=FA*CONT
	COSB=COS(X)
      SINB=SIN(X)
	COSB2=COSB*COSB
	COSB3=COSB2*COSB
	COSB4=COSB2*COSB2
C
	L2=L*L
	L3=L2*L
	L4=L2*L2
	ETA=EP*COSB
	ETA2=ETA*ETA
	T=TAN(X)
	T2=T*T
	N=A/SQRT(1.0-E2*SINB**2)
	X=SM+N*L2*SINB*COSB/2.0+N*L4*SINB*COSB3*(5.0-T2+9.0*ETA2
     *   +4.0*ETA2*ETA2)/24.0	
      Y=N*(COSB*L+(1.0-T2+ETA2)*COSB3*L3/6.0)
C	WRITE(*,*)' X,Y= ',X,Y
	RETURN
	END

      SUBROUTINE CURAY(JO)
	include "stack.inc"
      COMMON/STUFF/C(ICC),S(ICC),D(ICC),TH(ICC)
      COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC),NW,NWI(nwmax),
     &NWID(nwmax)
      COMMON/DEPTHS/DEPTH(ICC)
  8   DEPTH(1)=TTH(1)/2.0
      DO 10 J = 2,JO
 10    DEPTH(J)=DEPTH(J-1)+(TTH(J)+TTH(J-1))/2.
      DO 5 J = 1,JO
      Q = 6371.0 / (6371.0-DEPTH(J))
C      IF(FLAT) Q=1.
      CC(J)=(C(J))*Q
      SS(J)=(S(J))*Q
      DD(J)=(D(J))*Q
      TTH(J)=(TH(J))*Q
 5    CONTINUE
  9   RETURN
      END


            SUBROUTINE FIND (JO,PO,XO)
CL 07/26/10: correct for both free-surface reflection and discontinuity conversion
	include "stack.inc"
	DIMENSION XO(*)
      COMMON/TRAVEL/LP(ICC,2),LS(ICC,2) !07/26/10
      COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC),NW,NWI(nwmax),
     &NWID(nwmax)
      REAL*8 E(ICC),F(ICC),SXX,SYY,CXX,CYY,P
      K=JO
      P = PO
      PSQ = P ** 2
      DO 10 J = 1,K
	 SXX=SS(J)*P
	 SYY=CC(J)*P
	 CXX=DABS(1.-SXX*SXX)
	 CYY=DABS(1.-SYY*SYY)
      F(J) =SXX/DSQRT(CXX)
      E(J) =SYY/DSQRT(CYY)
 10   CONTINUE
      TOTEM = 0.0
 	DO 12 I=1,NW
	I1=1
	IF (I.GT.1) I1=1+NWID(I-1)
	NI=NWID(I)
	DO 11 J=I1,NI
c 11	TOTEM=TOTEM+E(J)*TTH(J)*LP(J)+F(J)*TTH(J)*LS(J)
 11     TOTEM=TOTEM+E(J)*TTH(J)*LP(J,2)+F(J)*TTH(J)*LS(J,2)  ! 07/26/10

12	XO(I)=TOTEM
      RETURN
      END


      SUBROUTINE FIND2 (JO,PO,TTO)
CL 07/26/10: correct for both free-surface reflection and discontinuity conversion
	include "stack.inc"
	DIMENSION TTO(*)
      COMMON/TRAVEL/LP(ICC,2),LS(ICC,2) !07/26/10
      COMMON/ORST/CC(ICC),SS(ICC),DD(ICC),TTH(ICC),NW,NWI(nwmax),
     &NWID(nwmax)
      REAL*8 E(ICC),F(ICC),PSQ,XX,YY,P
      K=JO
      P = PO
      PSQ = P ** 2
      DO 10 J = 1,K
	 XX=DABS(1./SS(J)-P)
	 YY=DABS(1./SS(J)+P)
      F(J) =DSQRT(XX*YY)
	 YY=DABS(1./CC(J)+P)
	 XX=DABS(1./CC(J)-P)
      E(J) =DSQRT(XX*YY)
 10   CONTINUE
      TOTEM = 0.0
 	DO 12 I=1,NW
	I1=1
	IF (I.GT.1) I1=1+NWID(I-1)
	NI=NWID(I)
	DO 11 J=I1,NI
c 11	TOTEM=TOTEM+E(J)*TTH(J)*LP(J)+F(J)*TTH(J)*LS(J)
 11     TOTEM=TOTEM+E(J)*TTH(J)*(LP(J,1)+LP(J,2))
     ^ +F(J)*TTH(J)*(LS(J,1)+LS(J,2))                        ! 07/26/10
12	TTO(I)=TOTEM
      RETURN
      END


C*************************************************************
       subroutine ichoose(i0,ikx,np,ip1,ip2,ik)
      
       integer ikx(np),i0
      
       npp=np-1
       if(i0.lt.ikx(ip1))then
         ik=ip1-1
       elseif(i0.ge.ikx(ip2))then
         ik=ip2
       else
         do ip=ip1+1,ip2
           if(i0.lt.ikx(ip)) then
	     ik=ip-1
	     return
	   endif
         enddo
       endif
       
       return
       end
C*************************************************************

C*************************************************************
       subroutine choose(x0,kx,np,ip1,ip2,ik)
      
       real kx(np),x0
      
       npp=np-1
       if(x0.lt.kx(ip1))then
         ik=ip1-1
       elseif(x0.ge.kx(ip2))then
         ik=ip2
       else
         do ip=ip1+1,ip2
           if(x0.lt.kx(ip)) then
	     ik=ip-1
	     return
	   endif
         enddo
       endif
       
       return
       end
C*************************************************************
