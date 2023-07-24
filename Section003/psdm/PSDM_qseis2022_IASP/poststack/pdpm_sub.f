

        SUBROUTINE GET_KX( nx, dkx, kx, kx2 )
        INTEGER  nx
        REAL     dkx
        REAL     kx( nx ), kx2( nx )

        nxhalf = nx / 2
        nx0    = nxhalf + 1
        nxp2   = nx + 2

        kx(1)=0.
        kx2(1)=0.

        do ix=2,nxhalf
          kx(ix)=(ix-1)*dkx
          kx2(ix)=kx(ix)*kx(ix)
        enddo

        kx(nx0)=-nxhalf*dkx
        kx2(nx0)=kx(nx0)*kx(nx0)

	do ix=2,nxhalf
	  kx(nxp2-ix)=-kx(ix)
	  kx2(nxp2-ix)=kx2(ix)
	enddo

 
        return
        end


        SUBROUTINE GET_OMEGA( nf, dw, w )
        INTEGER  nf
        REAL     dw
        REAL     w( nf )

        nfhalf = nf / 2
        nf0    = nfhalf + 1
        nfp2   = nf + 2

        do ifreq=1,nfhalf
          w(ifreq)=(ifreq-1)*dw
        enddo

        do ifreq=2,nfhalf
          w(nfp2-ifreq)=-w(ifreq)
        enddo

        w(nf0)=-nfhalf*dw

        return
        end

        
        SUBROUTINE CAL_KZ( nx, k02, kx2, kz )
        INTEGER  nx
        REAL     k02
        REAL     kx2( nx )
        REAL     kz( nx )
     
        nxhalf = nx / 2
        nx0    = nxhalf + 1
        nxp2   = nx + 2

	do ix=1,nx0
 	   if (k02.gt.kx2(ix)) then
 	      kz(ix)=sqrt(k02-kx2(ix))
 	   else
 	      kz(ix)=(0.,0.)
 	   endif
	enddo

	do ix=2,nxhalf
	   kz(nxp2-ix)=kz(ix)
	enddo

        return
        end

        
        SUBROUTINE CAL_PBS( nx, dz, kz, pbs )
        INTEGER  nx
        REAL     dz
        REAL     kz( nx )
        COMPLEX  im, cnum
        COMPLEX  pbs( nx ) 
 
        im = cmplx(0., 1.)
     
        nxhalf = nx / 2
        nx0    = nxhalf + 1
        nxp2   = nx + 2

        do ix=1,nx0
 	  if (kz(ix).gt.0) then
	    cnum = im * kz(ix) * dz
            pbs(ix) = cexp( cnum )
 	  else
 	    pbs(ix)=(0.,0.)
 	  endif
        enddo

	do ix=2,nxhalf
	   pbs(nxp2-ix)=pbs(ix)
	enddo

        return
        end


        SUBROUTINE SPATAP( nx, taperx, uzs )
        INTEGER  nx
        REAL     taperx( nx )
        COMPLEX  uzs( nx )

	  do ix=1,nx
	    uzs(ix) = uzs(ix)*taperx(ix)
        enddo

        return
        end


        SUBROUTINE CAL_REFV(nxx, nzz, nxs_vel, nxe_vel, nz,
     &                      irefvel, vel, v0)
        INTEGER  irefvel
        INTEGER  nxx, nzz
        INTEGER  nxs_vel, nxe_vel
        INTEGER  nz
        REAL     vel(nxx, nzz)
        REAL     v0(nz)

        nx = nxe_vel - nxs_vel + 1
C
C... minimum reference vel
        if (irefvel.eq.0) then
          do iz=1,nz
            vmin=200000.
            do ix=nxs_vel,nxe_vel
              if (vel(ix,iz).le.vmin) vmin=vel(ix,iz)
            enddo
            v0(iz)=vmin
          enddo
        endif

C
C... average reference vel
        if (irefvel.eq.1) then
          do iz=1, nz
            vsum=0.
            do ix=nxs_vel,nxe_vel
              vsum=vsum+vel(ix,iz)
            enddo
            v0(iz)=vsum/float(nx)
          enddo
        endif

C
C... min+ave reference vel
        if (irefvel.eq.2) then
          do iz=1, nz
            if (iz.le.300) then
              vmin=200000.
              do ix=nxs_vel,nxe_vel
                if (vel(ix,iz).le.vmin) vmin=vel(ix,iz)
              enddo
              v0(iz)=vmin
            else
              vsum=0.
              do ix=nxs_vel,nxe_vel
                vsum=vsum+vel(ix,iz)
              enddo
              v0(iz)=vsum/float(nx)
            endif
          enddo
        endif

        return
        end


        SUBROUTINE SRC_SPCTR(nf, wavelet, wlt)
        INTEGER nf, nnleng(1)
        REAL    wavelet( nf )
        COMPLEX wlt( nf ),wkspace(2*nf)
 
	do i=1, nf
	  wlt(i)=cmplx(wavelet(i),0.)
	enddo
C
C... FFT into frequency domain for the source 
*        call cwpfft(wlt, nf, -1)
        nnleng(1)=nf
	  wkspace(1)=(0.,0.)
        call fourt(wlt,nnleng,1,-1,1,wkspace)

        return
        end
 

        SUBROUTINE SRC_SHFT(nf, t0s, w, wlt)
        INTEGER  nf
        REAL     t0s
        REAL     w( nf )
        COMPLEX  im, cnum
        COMPLEX  wlt( nf )
 
        im = cmplx(0., 1.)

*        t0sr=t0s/1.5
        t0sr=t0s
        write(*,*)
        write(*,*)'t0s, t0sr=',t0s,t0sr

        if (t0sr.ne.0.) then
          do ifreq=1,nf
            omega= w(ifreq)
            cnum = im * omega * t0sr
*            wlt(ifreq)=-wlt(ifreq)* exp( cnum )
            wlt(ifreq)=wlt(ifreq)* cexp( cnum )
          enddo
        endif

        return
        end



	SUBROUTINE MUTE(ntmod, nxtemp, nxs_trin, ix, wavelet)
	INTEGER  ntmod
	INTEGER  nxtemp, nxs_trin
	INTEGER  ix
	REAL     wavelet(ntmod)
* mute
        ntm1=350
        ntm0=30
        dtmu=float(ntm1-ntm0)
        slope=dtmu/float(nxtemp)

        ixx=ix-nxs_trin+1
        if (ixx.lt.nxtemp) then
          itt=ntm1-ixx*slope
          do it=1,itt
	      wavelet(it)=0.
	    enddo
	  endif

        return
	  end


	SUBROUTINE GEOM_CORR(nf, dt, wavelet)
	INTEGER  nf
	REAL     wavelet(nf)

        do it=1,nf
          time=it*dt
          wavelet(it)=wavelet(it)*(time**0.8)
        enddo

        return
        end


        SUBROUTINE INIT_3DC(nx, nz, nf, a)
        INTEGER  nx, nz, nf
        COMPLEX  a(nx,nz,nf)

        do i=1,nf
          do iz=1,nz
	    do ix=1,nx
            a(ix,iz,i)=(0.,0.)
	    enddo
          enddo
        enddo

        return 
        end
        
	SUBROUTINE INIT_2DC(nx, nz, a)
        INTEGER  nx, nz
        COMPLEX  a(nx,nz)

        do ix=1,nx
          do iz=1,nz
            a(ix,iz)=(0.,0.)
          enddo
        enddo

        return 
        end


        SUBROUTINE INIT_1DR(nt, a)
        INTEGER  nt
        REAL     a( nt )

        do it=1,nt
         a(it)=0.
        enddo
         
        return
        end


