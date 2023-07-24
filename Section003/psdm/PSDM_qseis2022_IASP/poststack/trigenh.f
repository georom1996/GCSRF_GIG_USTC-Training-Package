
C-----------------------------------------------------------------------
C
C    trigenh:
C
C    computes tridiag matrix elements for left and right sides
C    of layer downward continuation equations given velocities
C    and basic geometry info.
C
C    a & b are associated with left side; aa & bb with right.
C
C
C     array dimension requirement:
C
C     complex
C         za -  nx
C         zb -  nx
C        zaa -  nx
C        zbb -  nx
C
C    complex scalars for absorbing boundaries
C
C      zbnx1 -  upper left  corner element, left  hand side
C      zbnxn -  lower right corner element, left  hand side
C     zbbnx1 -  upper left  corner element, right hand side
C     zbbnxn -  lower right corner element, right hand side
C
C
C
C  on entry:
C		wl   = angular frequency 		real, scalar
C		v    = velocity array for z-layer 	real, vector
C		dzm  = delta-z step size for layer	real, scalar
C		nx   = number of x-locations		int,  scalar
C 		dx   = delta-x spacing   		real, scalar
C       frac = fraction of velocity to use
C              frac = 1.0 for full velocity
C       iupd = up/downgoing flag, 1=up, -1=down
C		idip = 15, 45 or 65     	int,  scalar
C
C  on return:
C		complex difference operator is initialized for
C		use in seictri.
C
C        i.e.-  za(nx),zb(nx),zaa(nx),zbb(nx),
C	 	zbnx1, zbnxn, zbbnx1, zbbnxn
C		are returned for the layer.
C
C----------------------------------------------------------------------

	subroutine trigenh(  wl,v,vbar,za,zb,zaa,zbb,
     & 				zbnx1,zbnxn,zbbnx1,zbbnxn,
     &              dzm, nx, dx, frac, iupd, idip  )

C  Arguments
	integer nx, iupd, idip
	complex*8 za(nx),zb(nx),zaa(nx),zbb(nx)
	complex*8 zbnx1,zbnxn,zbbnx1,zbbnxn
	real*4 wl,v(nx),dzm,dx,frac
C  Local variables
	complex*8 zcc
	real*4 eta, a, b, beta, pi, dz, xupd, fact, dx4, den, alpha, wldv
	real*4 c3r, c2a, c2b, c1i, c2r, c2i, c, cc
	integer ix
C
C
	data beta /0.101/
	data pi   /3.14159265/

      twopi=pi+pi

        if (idip.eq.45) then
	  a=0.75
	  b=0.25
        else
	  a=0.5
	  b=0.
        endif
	
C  set sign of delta z and viscosity eta based on up/down flag
	xupd = iupd
	dz = dzm * xupd
	eta = pi * xupd
	if(idip.gt.45) eta = twopi * xupd
C
C  set velocity fraction
	if (frac.eq.0.) then
	    fact = 1.0
	else
	    fact = frac
	endif
C
C  constants
	dx4 = .25*dx
	den=wl*wl+eta*eta
	alpha=dx*dx/dz
	c3r=a-b

C* hybrid migration
	c2a=a/dz*wl/den
	c2b=a/dz*eta/den

C
C     calculate the tri-d elements
C
	do ix = 1 ,nx
C
            ssv= 1. - vbar/v(ix)
	    wldv = wl/v(ix)
	    c1i  = wldv*alpha
	    c2r  = c2b*(fact*v(ix))
	    c2i  = c2a*(fact*v(ix))+c1i*beta
            
 	    if (abs(ssv).le.1.0e-2) then
 	    rnza=0.
 	    rnzaa=0.
 	    else
	    rnza=c2r+.5*c3r*ssv
	    rnzaa=c2r-.5*c3r*ssv
 	    endif

	    za(ix)  = cmplx(rnza,c2i)
	    zaa(ix) = cmplx(rnzaa,c2i)
	    zb(ix)  = cmplx(0.,c1i)-2.*za(ix)
	    zbb(ix) = cmplx(0.,c1i)-2.*zaa(ix)
C
	enddo

C
C  absorbing boundary conditions ********
C
	  c         = xupd * wl/(fact*v(1))*dx4
	  cc        = c*c
	  zcc       = cmplx(1.-cc,2.*c)/(1.+cc)
	  zbnx1     = zb(1)+zcc*za(1)
	  zbbnx1    = zbb(1)+zcc*zaa(1)
C
 	  c         = xupd * wl/(fact*v(nx))*dx4
	  cc        = c*c
	  zcc       = cmplx(1.-cc,2.*c)/(1.+cc)
	  zbnxn     = zb(nx)+zcc*za(nx)
	  zbbnxn    = zbb(nx)+zcc*zaa(nx)

C
	return
	end
C
