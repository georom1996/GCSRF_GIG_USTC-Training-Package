PROGRAM signaltonoise
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!
!  Calculate a signal to noise ratio in a SAC file
!
!    1) define a time window around the signal of interest (s1-s2)
!       code simply finds the peak amplitude in this time window.
!
!    2) define a time window to characterize the noise level (n1-n2)
!       code calculates the RMS noise level in this window
!
!    3) SNR is defined here as peak amplitude/RMS noise level
!
! michael.thorne@utah.edu
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!
USE sac_i_o
IMPLICIT NONE
REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: amplitudes
REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: x !=====:added by Zhangzhou1996@gmail.com
REAL(KIND=4) :: s1, s2, n1, n2
REAL(KIND=4) :: time
REAL(KIND=4) :: maxamp, maxtime, mean_amp, DC, SLOPE, FLN, sum_y,sum_xy, sum_xx, sum_x, del, art, brt
!REAL(KIND=4) :: sum_y,sum_xy, sum_xx, sum_x, delta, a, b
REAL(KIND=4) :: summ, nsumm, SNR
INTEGER(KIND=4), DIMENSION(1) :: iptr
INTEGER(KIND=4)    :: NN, ios, J, Jz
CHARACTER(LEN=112) :: ifile
CHARACTER(LEN=12)   :: input


!    --  READ USER INPUT  --
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!
NN = IARGC()
IF (NN < 7) THEN
  write(*,'(a)') "usage:  sacsnr ifile -n n1 n2 -s s1 s2"
  write(*,'(a)') "        ifile: input SAC file"
  write(*,'(a)') "        -n n1 n2, noise time window"
  write(*,'(a)') "        -s s1 s2, signal time window"
  write(*,'(a)') "        e.g., sacsnr foo.sac -n 1300. 1400. -s 1450. 1460."
  STOP
ENDIF

! Read name of input sac file
CALL GETARG(1, ifile)
OPEN(UNIT=1,FILE=ifile,STATUS='OLD',IOSTAT=ios)
IF (ios > 0) THEN
  write(*,*) "ERROR - Input file: '", TRIM(adjustl(ifile)), "' does not exist ..."
  CLOSE(1)
  STOP
ENDIF
CLOSE(1) 

! Read first group of time window inputs
CALL GETARG(2,input)
IF (input(1:2) == '-s') THEN
  CALL GETARG(3,input)
  READ(input,*) s1
  CALL GETARG(4,input)
  READ(input,*) s2
ELSEIF (input(1:2) == '-n') THEN
  CALL GETARG(3,input)
  READ(input,*) n1
  CALL GETARG(4,input)
  READ(input,*) n2
ELSE
  write(*,*) "ERROR - Input argument: '", TRIM(adjustl(input)), "' is unknown ..."
  STOP
ENDIF

! Read second group of inputs
CALL GETARG(5,input)
IF (input(1:2) == '-s') THEN
  CALL GETARG(6,input)
  READ(input,*) s1
  CALL GETARG(7,input)
  READ(input,*) s2
ELSEIF (input(1:2) == '-n') THEN
  CALL GETARG(6,input)
  READ(input,*) n1
  CALL GETARG(7,input)
  READ(input,*) n2
ELSE
  write(*,*) "ERROR - Input argument: '", TRIM(adjustl(input)), "' is unknown ..."
  STOP
ENDIF

write(*,*) "Calculating SNR on file '", TRIM(adjustl(ifile)), "' ..."
write(*,*) "  Signal window: ", s1, "-", s2, "(sec)"
write(*,*) "   Noise window: ", n1, "-", n2, "(sec)"



!    --  READ INPUT SAC FILES  --
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!

!read ifile
CALL rbsac(ifile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,      &
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,  &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,amplitudes)

IF (nvhdr /= 6) THEN
  write(*,*) "ERROR - File: '", TRIM(adjustl(ifile)), "' appears to be of non-native &
  &byte-order or is not a SAC file."
  STOP
ENDIF

! check that windows are within sacfile
IF (s1 < b) THEN
  write(*,*) "ERROR:  s1 < begin time of trace.  Exiting..." 
  STOP
ELSEIF (n1 < b) THEN
  write(*,*) "ERROR:  n1 < begin time of trace.  Exiting..." 
  STOP
ENDIF

!    -- FIND MAX AMPLITUDE FOR PEAK --
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!
time = b
maxamp = -999.
maxtime = 0.


allocate(x(npts))
!===added by Zhangzhou1996@gmail.com
DO Jz=1,npts
        x(Jz) = dble(Jz) / dble(npts)
ENDDO
sum_x  = sum(x(1:npts))
sum_xx = sum(x(1:npts)*x(1:npts))
sum_y  = sum(amplitudes(1:npts))
sum_xy = sum(x(1:npts)*amplitudes(1:npts))
del = npts * sum_xx - sum_x * sum_x
art = ( npts* sum_xy - sum_y * sum_x ) / del
brt = ( sum_y * sum_xx - sum_x * sum_xy ) / del
! remove trend
do Jz=1, npts
       amplitudes(Jz) = amplitudes(Jz) - art*x(Jz)-brt
end do


!==added by Yangfan deng
mean_amp = sum(amplitudes)/npts
DO J=1,npts
	amplitudes(J) = amplitudes(J) - mean_amp
ENDDO
!===added by Yangfan deng

!      FLN = FLOAT(npts)
!      DC = 0.0
!      SLOPE = 0.0
!
!      DO J=1,npts
!        DC = DC + amplitudes(J)
!        SLOPE = SLOPE + amplitudes(J)*FLOAT(J)
!      ENDDO

!      DC = DC/FLN
!      SLOPE = 12.0*SLOPE/(FLN*(FLN*FLN-1.0)) - 6.0*DC/(FLN-1.0)
!      FLN = DC - 0.5*(FLN+1.0)*SLOPE
!      DO J=1,npts
!        amplitudes(J) = amplitudes(J) - DC
!        amplitudes(J) = amplitudes(J) - FLOAT(J)*SLOPE - FLN
!      ENDDO

!===

DO J=1,npts
  IF (time >= s1 .AND. time <= s2) THEN
    Amplitudes(J) = abs (amplitudes(J))  ! added
    IF (Amplitudes(J) > maxamp) THEN
      maxamp = Amplitudes(J)
      maxtime = time
    ENDIF
  ENDIF
  time = time + delta
ENDDO
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!

!  -- FIND RMS IN NOISE WINDOW --
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!
time = b
summ = 0.0
nsumm = 0.0
DO J=1,npts
  IF (time >= n1 .AND. time <= n2) THEN
    summ = summ + (amplitudes(J))**2
    nsumm = nsumm + 1.
  ENDIF
  time = time + delta
ENDDO
summ = sqrt((summ/nsumm))

SNR = maxamp/summ
!:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:=====:!

write(*,*) "Writing out to file 'snr.txt'"
write(*,*) "    output format:"
write(*,*) "    File name, time max signal, amplitude max signal, amplitude noise, SNR"
write(*,*) TRIM(ADJUSTL(ifile)), maxtime, maxamp, summ, SNR
OPEN(UNIT=2,FILE='snr.txt')
write(2,*) TRIM(ADJUSTL(ifile)), maxtime, maxamp, summ, SNR
CLOSE(2)


END PROGRAM signaltonoise
