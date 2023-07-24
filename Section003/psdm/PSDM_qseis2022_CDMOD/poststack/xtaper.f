        
	SUBROUTINE TAP_FUNC(taper, nx, nxleft, nxright, iflag)
        INTEGER  nx, nxleft, nxright
	INTEGER  iflag
	REAL     taper( nx )

************************************
*     iflag = 1 (one-sided taper)
*             2 (two-sided taper)
************************************

      pi = 3.1415926
      ntpleft = nxleft - 1
      ntpright = nxright - 1
      do kk=1,nx
	taper(kk)=1.0
      enddo

c     if (ntpleft.le.3.or.ntpright.le.3) then
c       return
c     endif

      do kk1=1,nx
        if (kk1.le.ntpleft.and.ntpleft.gt.3) then
   	  aa = 0.5-0.5*cos(pi*(kk1-1)/ntpleft)
 	  taper(kk1) = aa
        endif
        if (kk1.gt.nx-nxright.and.ntpright.gt.3) then
	  aa = 0.5-0.5*cos(pi*(nx-kk1)/ntpright)
	  taper(kk1) = aa
        endif
      enddo

      return
      end


************************************
*     Hanning taper subroutine
************************************
      subroutine staper_my(taper,nx,ntaper,iflag)
      dimension taper(nx)

************************************
*     iflag = -1 (left-sided taper)
*              1 (right-sided taper)
*              2 (two-sided taper)
************************************

      pi=3.1415926
      ntp=ntaper-1

      if (ntaper.le.1) then
        return
      endif

      if (iflag.eq.2) then
      do kk1=1,nx
      if (kk1.le.ntaper) then
 	aa=0.5-0.5*cos(pi*(kk1-1)/ntp)
 	taper(kk1)=aa
      endif
      if (kk1.gt.nx-ntaper) then
	aa=0.5-0.5*cos(pi*(nx-kk1)/ntp)
	taper(kk1)=aa
      endif
      enddo
      endif

      if (iflag.eq.-1) then
      do kk1=1,nx
      if (kk1.le.ntaper) then
	aa=0.5-0.5*cos(pi*(kk1-1)/ntp)
	taper(kk1)=aa
      endif
      enddo
      endif

      if (iflag.eq.1) then
      do kk2=1,nx
      if (kk2.gt.nx-ntaper) then
	aa=0.5-0.5*cos(pi*(nx-kk2)/ntp)
	taper(kk2)=aa
      endif
      enddo
      endif
      
      return
      end


