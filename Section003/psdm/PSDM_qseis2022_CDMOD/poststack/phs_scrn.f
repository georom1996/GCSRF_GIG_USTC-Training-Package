
	
	SUBROUTINE PHS_SHIFT(uz,pbs,dz,nx,k0,isign,k0dz)
        INTEGER isign
        REAL    k0, k0dz
        REAL    dz
	  REAL    rnum
        COMPLEX im, imk0dz
        COMPLEX uz( nx )
        COMPLEX pbs( nx )
	COMPLEX cnum,wkspace(2*nx)
	dimension nnleng(1)
	real scale


        im=cmplx(0.,1.)
        imk0dz = im * k0 * dz * isign
	  k0dz = k0 * dz * isign

C
C... FFT into Kx_domain
*        call cwpfft(uz, nx, -1)
        nnleng(1)=nx
	  wkspace(1)=(0.,0.)
        call fourt(uz,nnleng,1,-1,1,wkspace)
C
C... propagation in reference media
          if (isign.eq.1) then
            do ix=1,nx
	      uz(ix) = uz(ix) * pbs(ix)
            enddo
          else if (isign.eq.-1) then
	      do ix=1,nx
	      uz(ix)= uz(ix) * conjg(pbs(ix))
            enddo
          endif

C
C... IFFT into x_domain
*          call cwpfft(uz, nx, +1)
        nnleng(1)=nx
	  wkspace(1)=(0.,0.)
        call fourt(uz,nnleng,1,1,1,wkspace)
	  scale=1./float(nx)
	  do ix=1,nx
	    uz(ix) = uz(ix) * scale
	  enddo

	  return
	  end


	SUBROUTINE PHS_COR(uz,vbar,v,nx,k0dz)
        REAL    k0dz
        REAL    vbar, v( nx )
	  REAL    rnum
        COMPLEX uz( nx )
C
C... interaction within screen
	  do ix=1,nx
	    a= vbar / v(ix) - 1.
          rnum = k0dz * a
          uz(ix) = uz(ix) * cmplx(cos(rnum), sin(rnum))
        enddo
          
	  return
	  end

