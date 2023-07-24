
*########################################################*
*                                                        *
*	HDPMIG ---- Poststack Depth Migration with       *
*	            Pseudo-Screen  Propagator            *
*                   by hybrid implementation             *
*                                                        *
*########################################################*

*
C*    Parameters:
C*  
C*    nx - trace numbers
C*    dx - trace interval (CDP interval)
C*    nz - depth number
C*    dz - depth size (extrapolation step)
C*    nt - time samples of stacked section
C*    dt - sample rate 
C*    fmin - minimum frequency in calculation
C*    fmax - maximum frequency in Hz
C*    vmig - migration velocity
C*
C*----------------------------------------------------------------
C*    History:
C*
C*         Written by:
C*              01/10/98  --  Shengwen Jin
C*         Revised by:
C*              05/11/2004 -- Ling Chen
C*
C*         References:
C*              Jon Claerbout, 1985, Imaging the Earth's Interior
C*              Zaitian Ma, 1989,  Seismic Imaging Techniques
C*              Ru-Shan Wu, 1994,  Wide-angle ...
C*
C*---------------------------------------------------------------
C*

	program hdpmig

	parameter (nxx=2200,
     &	           nzz=1700,
     &	           ntt=4096)

	real w(ntt),wavelet(ntt),taperx(nxx)
	real k0,k02,k0dz

	complex us(nxx,ntt)
	complex uzs(nxx)
	complex work(ntt)
	complex pbs(nxx)
	complex image(nxx,nzz)

	character*100 inputdata
	character*100 modvelocity
	character*150 tx_data,image_data
	character*100 charatemp

	real timearray(2)

      real v(nxx),v0(nzz)
      real kx(nxx),kx2(nxx),kz(nxx)
      real vmig(nxx,nzz)

	complex za(nxx), zb(nxx), zaa(nxx), zbb(nxx), ze(nxx)
      complex zbnx1, zbnxn, zbbnx1, zbbnxn
                    
	integer intt(nxx),nnleng(1)
	character*100 trace_data
	real amp(ntt)
	complex wkspace(2*ntt)

	write(*,*)'----- Migration -----'
	write(*,*)
        inputdata='hdpmig.in'
	write(*,'(1x,"->",2x,a)')inputdata

      open(1,file=inputdata,status='old')

******************************************
*                                        *
*          input parameters              *
*                                        *
******************************************

*<2>
        write(*,*)
        write(*,*) 'imethod =?, iflag of reference vel =? '
C* imethod = 0: phase-shift; = 1: phase-screen; else: pseudo-screen
C* irefvel = 0: minimum velocity; else: average velocity 
        read(1,'(a)') charatemp
	read(1,*) imethod,irefvel, vscale
	write(*,*) '-> imethod =',imethod,' irefvel =',irefvel
	write(*,*) '-> vscale =',vscale

*<3>
	write(*,*)
	write(*,*)'Minimum and maximum frequencies=?'
	read(1,'(a)')charatemp
	read(1,*)fmin,fmax,ifreindl,ifreindr
	write(*,*)'->  Minimum frequency=',fmin
	write(*,*)'->  Maximum frequency=',fmax
	write(*,*)'-> ifreindl =',ifreindl,' ifreindr =',ifreindr
	ifreqint = 1
	ifreqint1 = ifreqint * 20

*<4,5>
c        write(*,*)'nxmod,nzmod,nx,nz=?'
	read(1,'(a)')charatemp
        read(1,*)nxmod,nzmod,nx,nz
	read(1,'(a)')charatemp
        read(1,*)dx,dz
	nx=max(nxmod,nx)
	nz=min(nzmod,nz)
	write(*,*)'->  nx=',nx,' nz=',nz,' dx=',dx,' dz=',dz
*<6>
c        write(*,*)'nr,nt,dt=?'
cl* nt0: fft length; ntb: samples of 1 - ntb-1 in each trace are muted
	read(1,'(a)')charatemp
      read(1,*) ntrace,nt,dt,nt0,ntb
      write(*,*)'->  ntrace=',ntrace,'  nt=',nt,' dt=',dt
	ntb=max(ntb,1)
	ntb=min(ntb,nt)

*<7>
      read(1,'(a)')charatemp
      read(1,*) idip
      write(*,*)'->  idip =', idip

	read(1,'(a)')charatemp
      read(1,*) nxleft, nxright
      write(*,*)'->  nxleft =', nxright

*<8>
	write(*,*)
	write(*,*)'ifmat = ?'
* ifmat: (=0: ascii, one file containing vp and vs, 1D model; 
*         else: binary, 2D model containing migration velocity)
	read(1,'(a)')charatemp
        read(1,*) ifmat
	write(*,*) '-> ifmat =',ifmat
	if(ifmat.eq.0) then
	  imethod = 0
	endif
*<9>
	write(*,*)
	write(*,*)'Model velocity data file=?'
	read(1,'(a)')charatemp
        read(1,'(a)') modvelocity
	write(*,'(1x,"->",2x,a)') modvelocity
*<10>
	write(*,*)
	write(*,*)'Reflected wave data file=?'
	read(1,'(a)')charatemp
        read(1,'(a)')tx_data
	write(*,'(1x,"->",2x,a)')tx_data
*<11>
	write(*,*)
	write(*,*)'Migrated image data file=?'
	read(1,'(a)')charatemp
        read(1,'(a)')image_data
	write(*,'(1x,"->",2x,a)')image_data
*<12>
c	write(*,*) 'selected trace interval intrace=?'
	read(1,'(a)') charatemp
	read(1,*) intrace
	if(intrace.eq.0)intrace=1
	write(*,*) 'intrace=',intrace
	if(intrace.lt.0)then
c	write(*,*)'trace index data file =?'
	  read(1,'(a)') charatemp
	  read(1,'(a)') trace_data
	else
c	write(*,*)'first trace index itrfirst = ?'
	  read(1,'(a)') charatemp
	  read(1,*) itrfirst
	endif

	intrace0=iabs(intrace)

	close(1)

* end of input parameters
* =======================

*+++++++++++++++++++++++++++++++++++*
*     Read velocities of the model  *
*+++++++++++++++++++++++++++++++++++*
C* backgound velocity indexes
	leftadd=(nx-nxmod)/2
	if(ifmat.ne.0) then
* read migration velocity from binary velocity data	  
	  nrecl=4*nzmod
          open(13,file=modvelocity,status='old',recl=nrecl,
     &	access='direct')

	    nxb=leftadd+1
	    nxe=leftadd+nxmod
          do ix=nxb,nxe
	      ix1=ix-nxb+1
	      read(13,rec=ix1) (vmig(ix,iz),iz=1,nzmod)
          enddo
          close(13)
C* extension
        do ix=1,leftadd
          do iz=1,nz
            vmig(ix,iz)=vmig(nxb,iz)
          enddo
        enddo
        do ix=nxe+1,nx
          do iz=1,nz
            vmig(ix,iz)=vmig(nxe,iz)
          enddo
        enddo
C
C... calculate background velocity within each slab
	  CALL CAL_REFV(nxx, nzz, 1, nx, nz, irefvel, vmig, v0)

	else
* 1D model: vp and vs are within one ascii file
        open(13,file=modvelocity,status='old')
  	  read(13,*)nlay
	  iz1=1
	  do 10 iz=1,nlay
	    read(13,*)vp,vs,des,tth,dep,ii
	    iz2=int(dep/dz+0.001)
	    iz2=min(iz2,nz)
	    if(iz2.ge.iz1)then
	      vv=(vp*vs)/(vp-vs)
	      do iz0=iz1,iz2
	        v0(iz0)=vv
	      enddo
	      iz1=iz2+1
	    endif
	    if(iz1.gt.nz)goto 11
10	  continue
11	  close(13)

	endif         ! end velocity input


C
C... recorded data
        nrecl=4*nt
        open(15,file=tx_data,status='unknown',recl=nrecl,
     &	access='direct')
C
C... Image result
	nrecl=4*nz
        open(16,file=image_data,status='unknown',recl=nrecl,
     &	access='direct')

*+++++++++++++++++++++++++++++++++++*
*      Parameter Initialization     *
*+++++++++++++++++++++++++++++++++++*

        frac=1.
	iupd=1
	pi=4.*atan(1.0)
	pid4=pi/4.
	tpi=pi+pi

C
C... get NF that can be FFT'ed,
*	call fftolen(ntmod, ntmod*2, nf)
* Ling
	nf = nt0
	nf0 = nf/2 + 1
	write(*,*) ' nf =', nf
C
C... calculate frequency smaple interval
	df = 1./(nf*dt)
	dw = tpi * df
C
C... calculate temporal frequency, w(nf)
	CALL GET_OMEGA( nf, dw, w )

C
C... get NX that can be FFT'ed,
*	call fftolen( nx, nx*2, nx)
*Ling: nx must be a power of 2
	write(*,*) ' nx =', nx

C
C... calculate wavenumber sample intervals
	ddkx = 1./(nx*dx)
	dkx  = tpi * ddkx

C
C... calculate horizontal wavenumber, kx(nx)
	CALL GET_KX( nx, dkx, kx, kx2 )

C
C... get taper function
        CALL TAP_FUNC( taperx, nx, nxleft, nxright, 2 )

C
C... calculate frequency index
	ifreqmax = min( nf0, ifix(fmax/df) + 1 )
	ifreqmin = max( 2  , ifix(fmin/df) + 1 )

	write(*,*)
	write(*,*)'ifreqmin,ifreqmax=',ifreqmin,ifreqmax,ifreqint
	write(*,*)'fmin,fmax=',(ifreqmin-1)*df,(ifreqmax-1)*df

	numfre=ifreqmax-ifreqmin+1
	if(ifreindl.gt.0.or.ifreindr.gt.0)then
 	  do it=1,numfre
	    amp(it)=1.
	  enddo
	  if(ifreindl.gt.0)then
	    ifreindl=min(ifreindl,numfre)
C*  freqency domain taper
          call staper_my(amp,numfre,ifreindl,-1)
	  endif
	  if(ifreindr.gt.0)then
	    ifreindr=min(ifreindr,numfre)
C*  freqency domain taper
          call staper_my(amp,numfre,ifreindr,1)
	  endif
	endif

	write(*,*)
	write(*,*)'fn=',1./(2.*dt)
	write(*,*)'df=',df,' dw=',dw,' dkx=',dkx

C
C... init image buffer
	CALL INIT_2DC(nxx, nz, image)
	CALL INIT_2DC(nxx, nf0, us)
C
C... init wavelet buffer
	CALL INIT_1DR(nt0, wavelet)
C
C... trace location
	if(intrace.lt.0)then
	  open(17,file=trace_data,status='old')
	  do ix=1,ntrace
	    read(17,*)intt(ix)
	  enddo
	  write(*,*)'intt(1:ntrace)='
	  write(*,*)(intt(ix),ix=1,ntrace)
	else

	  do ix=1,ntrace
	    ix1=(ix-1)*intrace0+itrfirst
	    intt(ix)=ix1

	  enddo

	endif

* end of initialization
* =====================

	cputime=dtime(timearray)
	write(*,*)'  -----CPU time=',cputime

*++++++++++++++++++++++++++++++++*
*     input reflected data       *
*++++++++++++++++++++++++++++++++*

	print*,'nx:',nx,' nt:',nt
        nnleng(1)=nt0
	wkspace(1)=(0.,0.)

	do ix=1,ntrace
	   ix1=intt(ix)
	   read(15,rec=ix)(wavelet(it),it=1,nt)
 
        
    
           do it=1,nt0
	       if(it.lt.ntb)then
	         work(it)=(0.,0.)
	       else
               work(it)=cmplx(wavelet(it),0.)
	       endif
           enddo



* fft over t
           call fourt(work,nnleng,1,-1,1,wkspace)


	     if(ifreindl.gt.0.or.ifreindr.gt.0)then
	       do it=1,numfre
	         ifreq=it+ifreqmin-1
	         work(ifreq)=work(ifreq)*amp(it)
	       enddo
	     endif

	   do ifreq=ifreqmin,ifreqmax

	      us(ix1+leftadd,ifreq)=work(ifreq)

	   enddo
	enddo

***************************************
*                                     *
*       Loop over frequency           *
*                                     *
***************************************

      write(*,*)' '
      write(*,*)'Starting the freqency loop ...'
	total=0.
	print*,'ifreqmin,ifreqmax,ifreqint:',ifreqmin,ifreqmax,ifreqint

	do 300 ifreq=ifreqmin,ifreqmax,ifreqint
 
        omega=w(ifreq)

*-----------------------------------*
*     Transform into Kx domain      *
*-----------------------------------*

	do ix=1,nx
	   uzs(ix)=us(ix,ifreq)
	enddo

**********************************
*                                *
*     downward propagation       *
*                                *
**********************************


	do 400 iz=1,nz

C
C... reference wavenumber
        vbar = v0( iz ) * vscale
        k0=omega/vbar
        k02=k0*k0

C
C... read migration velocity slice
	if(imethod.ne.0)then
        do ix=1,nx
          v(ix) = vmig(ix,iz)
        enddo
	endif

C
C... calculate vertical wavenumber, kz(nx)
	CALL CAL_KZ( nx, k02, kx2, kz )

C
C... calculate propagator, pbs(nx)
	CALL CAL_PBS( nx, dz, kz, pbs )

C
C... Spatial taper enforcement
	CALL SPATAP( nx, taperx, uzs )

C
C... Phase-shift Operator
	CALL PHS_SHIFT(uzs,pbs,dz,nx,k0,+1,k0dz)

C... Phase-Screen Correction
	if(imethod.ne.0)CALL PHS_COR(uzs,vbar,v,nx,k0dz)
	
      if (imethod.ne.0.and.imethod.ne.1) then
C
C... Wide-angle compensation term 
        CALL WIDE_CPNS(uzs, omega, v, vbar, dz, nx, dx, frac, iupd, 
     +                 za, zb, zaa, zbb, ze, idip)

      endif

*++++++++++++++++++++++++++++++++*
*     apply imaging condition    *
*++++++++++++++++++++++++++++++++*
        do ix=1,nx
          image(ix,iz)=image(ix,iz)+uzs(ix)
        enddo


********************************
*    Loop to next screen       *
********************************
400	continue

	if (mod(ifreq-ifreqmin,ifreqint1).eq.0) then
	   sum=0.
	   do iz=1,nz
	   do ix=1,nx
	      sum=sum+real(image(ix,iz))
	   enddo
	   enddo

	   write(*,"('  freq=',i5,'    sum(image)=',g11.4)")
     &	   ifreq,sum
	endif

        cputime=dtime(timearray)
	total=total+cputime
        write(*,*)'ifreq=',ifreq,ifreqmax,
     & ' --- CPU time, total CPU = ',cputime,total

********************************
*   Loop to next frequency     *
********************************
300	continue

* end of the freqency loop
* ========================
        write(*,*)'It is the end of the frequency loop.'
	write(*,*)
        write(*,*) '------- Total CPU time (sec.) =',total

***********************************
*                                 *
*     output migrated image       *
*                                 *
***********************************

* normalized image
	amax=0.
	do iz=1,nz
	do ix=1,nx
	   a=abs(real(image(ix,iz)))
	   if (a.gt.amax) amax=a
	enddo
	enddo
	write(*,*)'Max(real(image(ix,iz)))=',amax

	do ix=1,nxmod
	   ix1=ix+leftadd
	   write(16,rec=ix)(real(image(ix1,iz))/amax,iz=1,nz)
	enddo

        write(*,*)'Saving image in: ',image_data

	close (15)
	close (16)
	write(*,*)
	write(*,*)'output image size (nz,nx):',nz,' *',nxmod

	stop
	end
