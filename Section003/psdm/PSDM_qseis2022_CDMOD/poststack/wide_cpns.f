        SUBROUTINE WIDE_CPNS(uz, omega, v, vbar, dz, nx, dx, frac, iupd,
     +                       za, zb, zaa, zbb, ze, idip)			

        INTEGER nx
	INTEGER IUPD
	REAL    omega
	REAL    v(nx), vbar
	REAL    dx, dz
	COMPLEX uz(nx)
	COMPLEX zbnx1, zbnxn, zbbnx1, zbbnxn
	COMPLEX za(nx),zb(nx),zaa(nx),zbb(nx),ze(nx)

	call trigenh(omega,v,vbar,za,zb,zaa,zbb,
     +              zbnx1,zbnxn,zbbnx1,zbbnxn,
     +              dz, nx, dx, frac, iupd, idip)

        call trisolv(uz, za, zb, zaa, zbb,
     +              zbnx1, zbnxn, zbbnx1, zbbnxn,
     +              ze, nx)
        
	return
	end

