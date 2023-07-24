C	character kkdate*20
C	character kktime*20
C	ncdate = 20
C	call kadate(1998,133,ncdate,kkdate,nerr)
C	write(6,*)kkdate
C	call kadate(1998, 23,ncdate,kkdate,nerr)
C	write(6,*)kkdate
C	call kadate(1998,  3,ncdate,kkdate,nerr)
C	write(6,*)kkdate
C      call katime(23,33,43,890,12,kktime,nerr)
C	write(6,*)kktime
C      call katime( 3, 3, 3, 90,12,kktime,nerr)
C	write(6,*)kktime
C      call katime( 3, 3, 3,  0,12,kktime,nerr)
C	write(6,*)kktime
C      end

      subroutine katime(ihour,imin,isec,imsec,nctime,kktime,nerr)
*=====================================================================
* PURPOSE: To convert four integer fields representing hour, minute,
*          second, and millisecond to an alphanumeric equivalent
*          of the form:  HH:MM:SS.SSS
*=====================================================================
* INPUT ARGUMENTS:
*    IHOUR:   Integer hour field.
*    IMIN:    Integer minute field.
*    ISEC:    Integer second field.
*    IMSEC:   Integer millisecond field.
*    NCTIME:  Maximum length of output alphanumeric field.
*             NCTIME must be at least 12.
*=====================================================================
* OUTPUT ARGUMENTS:
*    KKTIME:  Equivalent alphanumeric field.
*             Set to all asterisks if an error occurs.
*    NERR:    Error flag. Set to 0 if no error occurred.
*             Potential error numbers: 0905, 0907.
*===================================================================
	integer ihour, imin, isec, imsec,nctime,nerr
	character kktime*(*)
	kktime(1:nctime) = ' '
	nerr = 0
	sec = isec + 0.001*imsec
	write(kktime,'(i2,1x,i2,1x,f6.3)')ihour, imin, sec
	kktime(3:3) = ':'
	kktime(6:6) = ':'
	do 1000 i=1,12
		if(kktime(i:i).eq.' ')kktime(i:i) = '0'
 1000	continue
	return
	end

      subroutine kadate(iyear,ijday,ncdate,kkdate,nerr)
*=====================================================================
* PURPOSE: To convert two integer fields representing year and
*          julian day to an alphanumeric equivalent
*          of the form:  MMM DD (JJJ), YYYY
*=====================================================================
* INPUT ARGUMENTS:
*    IYEAR:   Integer year field.
*    IJDAY:   Integer day field.
*    NCDATE:  Maximum length of output alphanumeric field.
*             NCDATE must be at least 18.
*=====================================================================
* OUTPUT ARGUMENTS:
*    KKDATE:  Equivalent alphanumeric date field.
*             Set to all asterisks if an error occurred.
*    NERR:    Error flag. Set to 0 if no error occurred.
*             Potential error numbers:
*=====================================================================
	integer iyear, ijday, imsec, nctime,nerr
	character kkdate*(*)
	character kmonth*4
      dimension kmonth(12)
      data kmonth/'JAN ','FEB ','MAR ','APR ','MAY ','JUN ',
     1            'JUL ','AUG ','SEP ','OCT ','NOV ','DEC '/


	nerr = 0
	call mnthdy(iyear,ijday,imonth,iday)
	kkdate(1:ncdate) = ' '
	kkdate(1:4)=kmonth(imonth)
	if(iday.gt.10)then
		write(kkdate(5:6),'(i2)')iday
	else
		kkdate(5:5) = '0'
		write(kkdate(6:6),'(i1)')iday
	endif
	kkdate(7:14)=' (   ), '
	if(ijday.gt.99)then
		write(kkdate(9:11),'(i3)')ijday
	else if(ijday.le.99 .and .ijday.ge.10)then
		kkdate(9:9) = '0'
		write(kkdate(10:11),'(i2)')ijday
	else
		kkdate(9:10)='00'
		write(kkdate(11:11),'(i1)')ijday
	endif
	
	write(kkdate(15:18),'(i4)')iyear

	return
	end

	subroutine etoh(epoch,date,str)
c-----
c	convert from epoch time to human time
c
c	epoch	- R*8 time in seconds relative to 0000 1 Jan 1970
c	date	- I*4 Julian date
c	str	- C*  printable string
c-----
	real*8 epoch
	integer*4 date
	character str*(*)
	integer*4 diy, doy, hour, minute, year, month
	integer*4 day
	real*4 second
	real*8 seclft
	

	str=' '
	doy = (epoch/86400.0d+00)
	seclft = dmod(epoch, 86400.0d+00)
	hour = 0
	minute = 0
	second = 0.00
c-----
c	compute hours minutes seconds
c-----
	if(seclft .ne. 0.00d+00)then
c-----
c				before 1970 subtract and add a day
c-----
		if(seclft .lt. 0.0d+00)then
			doy = doy - 1
			seclft = seclft + 86400.00d+00
		endif
		hour = (seclft/3600.00d+00)
		seclft = dmod(seclft,3600.0d+00)
		minute = seclft/60.0d+00
		second = dmod(seclft,60.0d+00)
	endif

	if(doy .ge. 0)then
		year = 1970
 1000		continue
			diy =  leapdy(year)
			if(doy .lt. diy)go to 2000
			doy = doy - diy
			year = year + 1
		go to 1000
	else
		year = 1969
 1100		continue
			diy =  leapdy(year)
			doy = doy + diy
			if( doy .gt. 0 ) go to 2000
			year = year - 1
		go to 1100
	endif
 2000	continue
	doy = doy + 1
	date = year*1000 + doy
	call mnthdy(year,doy,month,day)
	write(str,110) year,doy,year,month,day,hour,minute,second
  110	format(i4,'/',i3,' ',i4,'/',i2,'/',i2,' ',i2,':',i2,':',f6.3)
c-----
c	guarantee that there are no blanks in the string str
c-----
	do 2100 i=1,32
		if(str(i:i).eq.' ')str(i:i)='0'
 2100	continue
c-----
c	except between date and time
c-----
	str( 9: 9) = ' '
	str(20:20) = ' '
	return
	end

	function leapdy(yr)
	integer*4 yr
	logical t1, t2, t3
	t1 = mod(yr,4).ne.0
	t2 = mod(yr,100).ne.0
	t3 = mod(yr,400).ne.0
	if( .not.t1 .and. t2)then
		isleap = 1
		leapdy = 366
	elseif( .not.t3)then
		isleap = 1
		leapdy = 366
	else
		isleap = 0
		leapdy = 365
	endif
	return
	end

	subroutine mnthdy(year,doy,month,day)
	integer*4 year, doy, month, day
	integer*4 i, dim, leap
	integer*4 dmnth(12)
	data dmnth/31,28,31,30,31,30,31,31,30,31,30,31/
	if(leapdy(year).eq.366)then
		leap = 1
	else
		leap = 0
	endif
	day = doy
	do 100 i=1,12
		month = i
		dim = dmnth(i)
		if(leap.eq.1 .and. i.eq.2)dim = dim + 1
		if(day .le.dim)goto 1000
		day = day - dim 
  100	continue
 1000	continue
	return
	end


	subroutine htoe(year,month,day,hour,minute,second,epoch)
c-----
c	convert calendar date to epoch time since January 1, 1970
c-----
c	year	- I*4	year
c	month	- I*4	month
c	day	- I*4	day
c	hour	- I*4	hour
c	minute	- I*4	minute c	second	- I*4	second
c	second	- R*4	seconds
c	epoch	- R*8	time in seconds relative to 00:00 01 Jan 1970
c-----
	integer*4 year, month, day, hour, minute, date, diy
	real*4 second
	real*8 epoch, dtoepo
	integer*4 daymon(12)
	data daymon/0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
     1		304, 334/
	diy = daymon(month) + day
	if(leapdy(year).eq.366 .and. month.gt.2)diy=diy+1
	date = 1000*year + diy
c	write(6,*)'date=',date
	epoch = dtoepo(date) + hour * 3600.0d+00 + 
     1		minute * 60.0d+00 +dble(second)
	return
	end

c-----
c	convert julian date to epoch time
c-----
	function dtoepo(date)
	real*8 dtoepo
	integer*4 date, diy, cnt, days

	cnt = date / 1000
	days = 0
	if (cnt .gt. 1970)then
 1000		continue
		cnt = cnt -1
		if(cnt.lt.1970)go to 2000
			days = days + leapdy(cnt)
		go to 1000
	else if (cnt .lt. 1970)then
 1100		continue
		if(cnt.ge.1970)goto 2000
			days = days - leapdy(cnt)
			cnt = cnt + 1
		go to 1100
	endif
 2000	continue
	diy = (date -1) / 1000
	diy = (date -1 ) -  1000*diy
c	write(6,*)'days=',days,' diy=',diy
	dtoepo = (days + diy) * 86400.0d+00
	return
	end
		
		


	function lgstr(str)
c-----
c	function to find the length of a string
c	this will only be used with strings without
c	intervening spaces.
C	thus the first blank will indicate the end of the string
c-----
	character str*(*)
	integer lgstr
	n = len(str)
	lgstr = 1
	do 1000 i=n,1,-1
		lgstr = i
		if(str(i:i).ne.' ' ) go to 100
 1000	continue
  100	continue
	return
	end
c-----
c	handler to catch floating point errros for SUN
c	do not use in production work, since it uses a
c	lot of cycles
c-----
c	do a f77 -g compile, and
c	use dbx, catch FPE, run
c	the dbx  will stop at floating point errors
c-----
      subroutine mchdep()
      return
      end

      subroutine handler(i1,i2,i3,i4)
      return
      end
	subroutine mgtarg(i,name)
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MGTARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c	return the i'th command line argument
c
c	This version works on SUN, IBM RS6000
c-----
	character name*(*)
		call getarg(i,name)
	return
	end
	function mnmarg()
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MNMARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
	integer mnmarg
		mnmarg = iargc()
	return
	end
	subroutine rsac1(infile,y,npts,btime,dt,maxpts,nerr)
	character infile*(*)
	real y(maxpts)
	integer npts
	real btime
	real dt
	integer maxpts
	integer nerr
c-----
c	PURPOSE READ AN EVENLY SPACED SAC FILE
c
c	read a binary sac file with evenly sampled data
c	infile	Char	name of file
c	y	R	array of values
c	npts	I	number of points in data
c	btime	R	start time
c	dt	R	sample interval
c	maxpts	I	maximum number of points to read
c	nerr	I	error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)

	npts = maxpts
c-----
c	rad up to maxpts, indicate the number actually read
c-----
	call brsac(1,maxpts,infile,y,nerr)
	call getnhv('NPTS    ',npts,ierr)
	call getfhv('DELTA   ',dt  ,ierr)
	return
	end

	subroutine wsac0(ofile,x,y,nerr)
	character ofile*(*)
	real x(*)
	real y(*)
	integer nerr

	logical leven
c-----
c	WRITE  SAC FILE USING CURRENT HEADER VALUES
c	however we should look at the header value LEVEN
c	to decide if we should write out at x,y pairs or as
c	y values only
c
c	RBH 2000 08 31 look at the header value,
c	then do a wsac1 or a wsac2
c
c	write a binary sac file with evenly sampled data
c
c	ofile	Char	name of file
c	ydummy	R	array of values
c	y	R	array of values
c	nerr	I	error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)

	npts = ihdr(10)
c-----
c	this is wrong
c	'To write a SAC file to disk using current header values'
C	well I am forcing even spaced data
c		had to do this to kludge owens SAC
c-----
	call getlhv('LEVEN   ',leven,nerr)
	call getnhv('NPTS    ',npts,nerr)
	call getfhv('DELTA   ',dt,nerr)
	call getfhv('B       ',b ,nerr)
	if(leven)then
		call wsac1(ofile,y,npts,b,dt,nerr)
	else
		call wsac2(ofile,x,npts,y,nerr)
	endif
C	call bwsac(1,npts,ofile,y)
	nerr = 0
	return
	end

	subroutine wsac1(ofile,y,npts,btime,dt,nerr)
	character ofile*(*)
	real y(npts)
	integer npts
	real btime
	real dt
	integer nerr 
c----- 
c	PURPOSE: WRITE AN EVENLY SPACED SAC FILE
c
c	write a binary sac file with evenly sampled data
c	ofile	Char	name of file
c	y	R	array of values
c	npts	I	number of points in data
c	btime	R	start time
c	dt	R	sample interval
c	maxpts	I	maximum number of points to read
c	nerr	I	error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c	create a new header
c-----
C	call newhdr()
c-----
c	place essential values into the new header
c-----
	call setnhv('NPTS    ',npts,nerr)
	call setfhv('DELTA   ',dt  ,nerr)
	call setfhv('B       ',btime  ,nerr)
	call setihv('IFTYPE  ','ITIME   ',nerr)
	e = btime + (npts -1 )*dt
	call setfhv('E       ',e     ,nerr)
	call setlhv('LEVEN   ',.true.,nerr)
	call setlhv('LOVROK  ',.true.,nerr)
	call setlhv('LCALDA  ',.true.,nerr)
	call bwsac(1,npts,ofile,y)
	nerr = 0
	return
	end

	subroutine wsac2(ofile,x,npts,y,nerr)
	character ofile*(*)
	real y(npts), x(npts)
	integer npts
	integer nerr 
c----- 
c	PURPOSE: WRITE AN UNEVENLY SPACED SAC FILE
c
c	write a binary sac file with evenly sampled data
c	ofile	Char	name of file
c	y	R	array of values
c	npts	I	number of points in data
c	x	R	array of independent variable
c	nerr	I	error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c	create a new header
c-----
	call newhdr()
c-----
c	place essential values into the new header
c-----
	call setnhv('NPTS    ',npts,nerr)
	call setihv('IFTYPE  ','ITIME   ',nerr)
	call setfhv('B       ',y(1)  ,nerr)
	call setfhv('E       ',y(npts)  ,nerr)
	call bwsac2(1,npts,ofile,x,y,npts)
	nerr = 0
	return
	end

	subroutine rsac2(ofile,x,npts,y,maxpts,nerr)
	character ofile*(*)
	real y(npts), x(npts)
	integer npts
	integer nerr 
	integer maxpts
c----- 
c	PURPOSE: WRITE AN UNEVENLY SPACED SAC FILE
c
c	write a binary sac file with evenly sampled data
c	ofile	Char	name of file
c	y	R	array of values
c	npts	I	number of points in data
c	x	R	array of independent variable
c	nerr	I	error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c	create a new header
c-----
	
c-----
c	place essential values into the new header
c-----
	call brsac2(1,maxpts,ofile,x,y,npts)
	nerr = 0
	return
	end

	subroutine brsac (IRU,LN,name,data,nerr)
c-----
c	IRU	I*4	logical unit for IO
c	LN	I*4	length of data array
c	name	C*	Name of file to be opened
c	rhdr	R*4	Real header
c	ihdr	I*4	Integer Header
c	chdr	C*	Character Header
c	data	R*4	Array of trace values
c	nerr	I*4	-1 file does not exist
c			-2 data points in file exceed dimension
c
c	NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c	LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----
	real*4 data(LN)

	logical ext

	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*(*) name
c-----
c  Read real and integer header blocks to find actual number
c  of waveform data points which is stored in ihdr(10).
c
c-----
	inquire(file=name,exist=ext)
	if(.not. ext)then
		ihdr(10) = 0
		nerr = -1
		return
	endif
		nerr = 0
		open (IRU,file=name,form='unformatted',
     &			access='direct',recl=440,status='old')
			read (IRU,rec=1) (rhdr(i), i=1,70),
     &					 (ihdr(i), i=1,40)
		close (IRU)
c-----
c
c  Read header and waveform data blocks using recored length of 158*4=632.
c
c-----
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
			nerr = -2
		else 
			maxpts = ihdr(10)
			nerr = 0
		endif
		nrec=632+4*maxpts
		nread = 0
c-----
C-----
C	The true simple code is commented out. The code to work
C	with an early buggy SUNOS 4.1 FORTRAN follows
C-----
C		nrec=632+4*ihdr(10)
C		open (IRU,file=name,form='unformatted',
C     &			access='direct',recl=nrec)
C			read (IRU,rec=1) (rhdr(i), i=1,70),
C     &					 (ihdr(i), i=1,40),
C     &					 (chdr(i), i=1,24),
C     &					 (data(i), i=1,ihdr(10))
C		close (IRU)
C-----
		nrec=632+4*maxpts
		nread = 0
c-----
c	because of SUNOS Fortran problems with IO transfers 
c	more than 2048 bytes, read these  chunks in 
c----- 
		ndat = maxpts
		if(nrec.gt.2048)then
			open (IRU,file=name,form='unformatted',
     &				access='direct',recl=2048)
			ndat1 = (2048 - 632) / 4
			irec = 1
			read (IRU,rec=irec,err=1001) (rhdr(i), i=1,70),
     &					 (ihdr(i), i=1,40),
     &					 (chdr(i), i=1,24),
     &					 (data(i), i=1,ndat1)
			nread = nread + ndat1
 1000			continue
			nl = nread + 1
			nh = nl + 512 - 1
			if(nh.gt.ndat)then
				nh = ndat
			endif
			if(nl.gt.ndat)go to 1001
			irec = irec + 1
			read (IRU,rec=irec,err=1001) (data(i), i=nl,nh)
			nread = nread + (nh-nl+1)

			go to 1000
 1001			continue
		close (IRU)
		else
			open (IRU,file=name,form='unformatted',
     &				access='direct',recl=nrec)
			read (IRU,rec=1) (rhdr(i), i=1,70),
     &					 (ihdr(i), i=1,40),
     &					 (chdr(i), i=1,24),
     &					 (data(i), i=1,ndat)
		close (IRU)
		endif
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
		else 
			maxpts = ihdr(10)
		endif
	return
	end

	subroutine arsac (IRU,LN,name,data,nerr)
c-----
c	IRU	I*4	logical unit for IO
c	LN	I*4	length of data array
c	name	C*	Name of file to be opened
c	rhdr	R*4	Real header
c	ihdr	I*4	Integer Header
c	chdr	C*	Character Header
c	data	R*4	Array of trace values
c	nerr	I*4	-1 file does not exist
c			-2 data points in file exceed dimension
c
c	NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c	LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  This routine reads files written in SAC alpha (ASCII) format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----
	logical ext
	real*4 data(LN)
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*(*) name
c-----
	inquire(file=name,exist=ext)
	if(.not. ext)then
		ihdr(10) = 0
		nerr = -1
		return
	endif
	nerr = 0
		open (IRU,file=name,status='old',access='sequential')
c----		rewind IRU
c  Read real header block.
c-----
		j1=1
		j2=5
		do 1110  i=1,14
			read (IRU,'(5g15.0)') (rhdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1110		continue
c-----
c  Read integer header block.
c-----
		j1=1
		j2=5
		do 1120 i=1,8
			read (IRU,'(5i10)') (ihdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1120		continue
c-----
c  Read character header block.
c-----
		j1=1
		j2=3
		do 1130 i=1,8
			read (IRU,'(3a8)') (chdr(j), j=j1,j2)
			j1=j1+3
			j2=j2+3
 1130		continue
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
			nerr = -2
		else 
			maxpts = ihdr(10)
			nerr = 0
		endif
c-----
c  Read waveform data organized in a five columns block.
c-----
		nrow=(maxpts/5)+1
		do 1140 i=1,nrow
			l=i*5
			k=l-4
			read (IRU,'(5g15.0)') (data(j), j=k,l)
 1140		continue
		close (IRU)
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
		else 
			maxpts = ihdr(10)
		endif
	return
	end

	subroutine getfhv(strcmd,fval,nerr)
c-----
c	Get float header value
c
c	strcmd	C*8	String to key on
c	val	R*4	Real value returned
c	nerr	I*4	Error condition
c				0 no error
c				1336 Value not defined
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output real header
c-----
	fval = -12345.
	nerr = -1
	do 1000 i=1,70
		if(streql(strcmd,rstr(i))) then
			nerr = 0
			fval = rhdr(i)
		endif
 1000	continue
	return
	end

	subroutine getnhv(strcmd,ival,nerr)
c-----
c	Get integer header value
c
c	str	C*8	String to key on
c	ival	R*4	integer value returned
c	nerr	I*4	Error condition
c				0 no error
c				1336 Value not defined
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output integer header
c-----
	ival = -12345
	nerr = -1
	do 2000 i=1,40
		if(streql(strcmd,istr(i))) then
			nerr = 0
			ival = ihdr(i)
		endif
 2000	continue
	return
	end

	subroutine getkhv(strcmd,cval,nerr)
c-----
c	Get character header value
c
c	strcmd	C*8	String to key on
c	cval	C*8	character value returned
c	nerr	I*4	Error condition
c				0  no error
c				1336 Value not defined
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql

	character cval*8
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output character header
c-----
	nerr = -1
	do 3000 i=1,24
		if(streql(strcmd,cstr(i))) then
			nerr = 0
			cval = chdr(i)
		endif
 3000	continue
	return
	end

	subroutine getlhv(strcmd,lval,nerr)
	character strcmd*(*)
	logical lval
	
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	input logical header
c-----
	call getnhv(strcmd,ival,nerr)
	if(ival.eq.0)then
		lval = .false.
	else
		lval = .true.
	endif
	return
	end

c---------------------------------------------------------
	subroutine bwsac (IWU,LN,name,data)
c---------------------------------------------------------

c
c  This routine writes out a waveform data in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character name*(*)
	real*4 data(LN)
c
c  The actual number of waveform data points is stored in integer
c  header 10. The file recored length is 158*4=632.
c
		nrec=632+4*ihdr(10)
		open (IWU,file=name,form='unformatted',
     &			access='direct',recl=nrec,status='replace')
c  Zhangzhou change the open status from replace into direct
c  Zhangzhou commoned the rewind IWU line
c		rewind IWU
		write (IWU,rec=1) (rhdr(i),i=1,70),
     &				  (ihdr(k),k=1,40),
     &				  (chdr(j),j=1,24),
     &				  (data(l), l=1,ihdr(10))
		close (IWU)
	return
	end
c---------------------------------------------------------
	subroutine awsac (IWU,LN,name,data)
c---------------------------------------------------------
c
c  This routine writes out files in SAC alpha (ASCII) format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
	real*4 data(LN)
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character name*(*)
c
		open (IWU,file=name,status='unknown',
     &				 access='sequential')
		rewind IWU
c
c  Write real header block.
c
		j1=1
		j2=5
		do 1100 i=1,14
			write (IWU,'(5g15.7)') (rhdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1100		continue
c
c  Write integer header block.
c
		j1=1
		j2=5
		do 1110 i=1,8
			write (IWU,'(5i10)') (ihdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1110		continue
c
c  Write character header block.
c
		j1=1
		j2=3
		do 1120 i=1,8
			write (IWU,'(3a8)') (chdr(j), j=j1,j2)
			j1=j1+3
			j2=j2+3
 1120		continue
c
c  Ensure the last row is padded with zeros, if actual number of
c  waveform points is less than the product of number of rows
c  and number of columns which constitutes the data block.
c
		nrow=(ihdr(10)/5)+1
		nrc5=nrow*5
		if (nrc5 .gt. ihdr(10)) then
			nrcx=ihdr(10)+1
			do 1140 i=nrcx,nrc5
				data(i)=0.0
 1140			continue
		end if
c
c  Write waveform data in five columns format.
c
		do 1150 i=1,nrow
			k=i*5
			j=k-4
			write (IWU,'(5g15.7)') (data(l), l=j,k)
 1150		continue
		close (IWU)
	return
	end

	subroutine setfhv(strcmd,fval,nerr)
c-----
c	Set float header value
c
c	strcmd	C*8	String to key on
c	fval	C*8	real value set
c	nerr	I*4	Error condition
c				0  no error
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output real header
c-----
	do 1000 i=1,70
		if(streql(strcmd,rstr(i))) rhdr(i) = fval
 1000	continue
	return
	end

	subroutine setnhv(strcmd,ival,nerr)
c-----
c	Set integer header value
c
c	strcmd	C*8	String to key on
c	ival	C*8	integer value set
c	nerr	I*4	Error condition
c				0  no error
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	integer ival, nerr

	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output integer header
c-----
	do 2000 i=1,40
		if(streql(strcmd,istr(i))) ihdr(i) = ival
 2000	continue
	return
	end

	subroutine setkhv(strcmd,cval,nerr)
c-----
c	Set character header value
c
c	strcmd	C*8	String to key on
c	cval	C*8	character value set
c	nerr	I*4	Error condition
c				0  no error
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	character cval*8
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output character header
c-----
	do 3000 i=1,24
		if(streql(strcmd,cstr(i))) chdr(i) = cval
 3000	continue
	return
	end


	subroutine setlhv(strcmd,lval,nerr)
	character strcmd*(*)
	logical lval
	
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output logical header
c-----
	if(lval)then
		call setnhv(strcmd,1,nerr)
	else
		call setnhv(strcmd,0,nerr)
	endif
	return
	end


	subroutine newhdr()
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	call inihdr()
	return
	end

	subroutine inihdr()
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	do 1100 i=1,70
		rhdr(i)= -12345.0
 1100	continue
	do 1110 i=1,35
		ihdr(i)= -12345
 1110	continue
	ihdr(7)=6
	ihdr(8)=0
	ihdr(9)=0
	do 1120 i=36,40
		ihdr(i)=0
 1120	continue
	do 1130 i=1,24
		chdr(i)='-12345  '
 1130	continue
	return
	end


	logical function streql(str1,str2)
	character str1*(*), str2*(*)
	character nstr1*8, nstr2*8
c-----
c	determine if two strings are equal
c-----
	nstr1 = ' '
	l1 = lgstr(str1)
	nstr1(1:l1) = str1(1:l1)
	nstr2 = ' '
	l2 = lgstr(str2)
	nstr2(1:l2) = str2(1:l2)
c-----
	if(nstr1 .eq. nstr2)then
		streql = .true.
	else
		streql = .false.
	endif
	return 
	end

	subroutine getihv(strcmd,strval,nerr)
c-----
c	Get enumerated header value
c
c	strcmd	C*8	String to key on
c	strval	C*8	real value set
c	nerr	I*4	Error condition
c			0  no error
c			1336 Header variable undefined
c			1337 Header variable does not exist
c-----
	character strcmd*(*), strval*8
	parameter (NEVAL=50)
	character*8 eval(NEVAL)
c-----
c	header integer equivalents of enumerated values
c	e.g., IDISP == 2
c-----
      data eval/'ITIME   ','IRLIM   ','IAMPH   ','IXY     ','IUNKN   ', 
     1	'IDISP   ', 'IVEL    ', 'IACC    ', 'IB      ', 'IDAY    ', 
     2	'IO      ', 'IA      ', 'IT0     ', 'IT1     ', 'IT2     ', 
     3	'IT3     ', 'IT4     ', 'IT5     ', 'IT6     ', 'IT7     ', 
     4	'IT8     ', 'IT9     ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5	'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6	'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7	'ISRO    ', 'INUCL   ', 'IPREN   ', 'IPOSTN  ', 'IQUAKE  ', 
     8	'IPREQ   ', 'IPOSTQ  ', 'ICHEM   ', 'IOTHER  ', 'IGOOD   ', 
     9	'IGLCH   ', 'IDROP   ', 'ILOWSN  ', 'IRLDTA  ', 'IVOLTS  '/
c-----
		call getnhv(strcmd,nval,nerr)
		strval = ' '
		if(nerr.eq.0)then
			if(nval.ge.1 .and. nval.le.NEVAL)then
				strval = eval(nval)
			endif
		endif
	return
	end

	subroutine setihv(strcmd,strval,nerr)
c-----
c	Set enumerated header value
c
c	strcmd	C*8	String to key on
c	strval	C*8	real value set
c	nerr	I*4	Error condition
c			0  no error
c			1336 Header variable undefined
c			1337 Header variable does not exist
c-----
	character strcmd*(*), strval*8
	parameter (NEVAL=50)
	character*8 eval(NEVAL)
	character*8 ival(NEVAL)
	logical streql
c-----
c	header integer equivalents of enumerated values
c	e.g., IDISP == 2
c-----
      data eval/'ITIME   ','IRLIM   ','IAMPH   ','IXY     ','IUNKN   ', 
     1	'IDISP   ', 'IVEL    ', 'IACC    ', 'IB      ', 'IDAY    ', 
     2	'IO      ', 'IA      ', 'IT0     ', 'IT1     ', 'IT2     ', 
     3	'IT3     ', 'IT4     ', 'IT5     ', 'IT6     ', 'IT7     ', 
     4	'IT8     ', 'IT9     ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5	'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6	'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7	'ISRO    ', 'INUCL   ', 'IPREN   ', 'IPOSTN  ', 'IQUAKE  ', 
     8	'IPREQ   ', 'IPOSTQ  ', 'ICHEM   ', 'IOTHER  ', 'IGOOD   ', 
     9	'IGLCH   ', 'IDROP   ', 'ILOWSN  ', 'IRLDTA  ', 'IVOLTS  '/
c-----
c	equivalence of header field position and enumerated value
c	e.g., IDEP can be IUNKN, IDISP, IVEL, IVOLTS or IACC
c-----
      data ival/'IFTYPE  ','IFTYPE  ','IFTYPE  ','IFTYPE  ','IDEP    ', 
     1	'IDEP    ', 'IDEP    ', 'IDEP    ', 'IZTYPE  ', 'IZTYPE  ', 
     2	'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 
     3	'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 
     4	'IZTYPE  ', 'IZTYPE  ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5	'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6	'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7	'ISRO    ', 'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 
     8	'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 'IQUAL   ', 'IQUAL   ', 
     9	'IQUAL   ', 'IQUAL   ', 'IQUAL   ', 'ISYNTH  ', 'IDEP    '/

c-----
c	now do the work, parse the table for the match
c		strcmd = ival and strval = eval, then 
c		using the table index, I, 
c			do a call setnhv(strcmd,I,nerr)
c
c	However, the IUNKN is used in both IDEP and IZTYPE 
c	and IOTHER is used in both IEVTYP and IQUAL
c	
c-----
	nerr = 0
	if(streql(strcmd,'IDEP    ')
     1			.and. streql(strval,'IUNKN   '))then
		call setnhv('IDEP    ',5,nerr)
	else if(streql(strcmd,'IZTYPE  ')
     1			.and. streql(strval,'IUNKN   '))then
		call setnhv('IZTYPE  ',5,nerr)
	else if(streql(strcmd,'IEVTYP  ')
     1			.and. streql(strval,'IOTHER  '))then
		call setnhv('IEVTYP  ',44,nerr)
	else if(streql(strcmd,'IQUAL   ')
     1			.and. streql(strval,'IOTHER  '))then
		call setnhv('IQUAL   ',44,nerr)
	else
		nerr = 1336
c-----
c	IFTYPE
c-----
		do 100 i=1,NEVAL
			if(
     1				streql(strval,eval(i)))then
				call setnhv(strcmd,i,nerr)
			endif
  100	continue
	endif
	return
	end
		
		
	
c---------------------------------------------------------
	subroutine bwsac2 (IWU,LN,name,x,y,npts)
c---------------------------------------------------------

c
c  This routine writes out a waveform data in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)

	character name*(*)
	real*4 x(LN)
	real*4 y(LN)
c
c  The actual number of waveform data points is stored in integer
c  header 10. The file recored length is 158*4=632.
c
		nrec=632+8*npts
		open (IWU,file=name,form='unformatted',
     &			access='direct',recl=nrec,status='unknown')
		rewind IWU
		write (IWU,rec=1) (rhdr(i),i=1,70),
     &				  (ihdr(k),k=1,40),
     &				  (chdr(j),j=1,24),
     &				  (x(l), l=1,npts),
     &				  (y(l), l=1,npts)
		close (IWU)
	return
	end
c---------------------------------------------------------
	subroutine brsac2 (IRU,LN,name,x,y,npts)
c---------------------------------------------------------
c-----
c	IRU	I*4	logical unit for IO
c	LN	I*4	length of data array
c	name	C*	Name of file to be opened
c	rhdr	R*4	Real header
c	ihdr	I*4	Integer Header
c	chdr	C*	Character Header
c	x	R*4	Array of x values
c	y	R*4	Array of y values
c	nerr	I*4	-1 file does not exist
c			-2 data points in file exceed dimension
c
c	NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c	LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----

	logical ext

	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)

	character*(*) name
	integer LN
	real*4 x(LN), y(LN)
	integer npts

	
c-----
c  Read real and integer header blocks to find actual number
c  of waveform data points which is stored in ihdr(10).
c
c-----
	inquire(file=name,exist=ext)
	if(.not. ext)then
		ihdr(10) = 0
		nerr = -1
		return
	endif
		nerr = 0
		open (IRU,file=name,form='unformatted',
     &			access='direct',recl=440,status='old')
			read (IRU,rec=1) (rhdr(i), i=1,70),
     &					 (ihdr(i), i=1,40)
		close (IRU)
c-----
c
c  Read header and waveform data blocks using recored length of 158*4=632.
c
c-----
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
			nerr = -2
		else 
			maxpts = ihdr(10)
			nerr = 0
		endif
		npts = maxpts
		nrec=632+8*npts
		open (IRU,file=name,form='unformatted',
     &			access='direct',recl=nrec,status='unknown')
		rewind IRU
		read (IRU,rec=1) (rhdr(i),i=1,70),
     &				  (ihdr(k),k=1,40),
     &				  (chdr(j),j=1,24),
     &				  (x(l), l=1,npts),
     &				  (y(l), l=1,npts)
	return
	end


	subroutine gcmdln(sacbin,fname)
	logical sacbin
	character fname*(*)

	character name*80
	nmarg = mnmarg()
	sacbin = .true.
	fname = ' '
	i = 0
   10	continue
		i = i + 1
		if(i.gt.nmarg)go to 11
		call mgtarg(i,name)
		if(name(1:2).eq.'-A')then
			sacbin = .false.
		else if(name(1:2).eq.'-B')then
			sacbin = .true.
		else if(name(1:2).eq.'-?')then
			call usage()
		else 
			fname = name
		endif
	go to 10
   11	continue
	return
	end

	subroutine usage()
	parameter (LER=6)
	write(LER,*)'Usage: shwsac [-B] [-A] [-?] file_name'
	write(LER,*)
     1		' -B	file is a binary SAC file (default)'
	write(LER,*)
     1		' -A	file is an ASCII SAC file (default)'
	write(LER,*)
     1		' -?	Command help - do not run'
	write(LER,*)
     1		' file_name Name of SAC file (required)'
	stop
	end

	subroutine parsesac(LN,x)

	parameter(LIN=5, LOT=6, LER=0)

	real*4 rval
	integer*4 ival
	character*8 cval
	logical lval

	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	real*4 x(LN)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA'  , 'DEPMIN'  , 'DEPMAX'   , 'SCALE'   , 'ODELTA'   , 
     1	'B'       , 'E'      , 'O'        , 'A'       , 'FMT'      , 
     1	'T0'      , 'T1'      , 'T2'      , 'T3'      , 'T4'      , 
     1	'T5'      , 'T6'      , 'T7'      , 'T8'      , 'T9'      , 
     1	'F'       , 'RESP0'   , 'RESP1'   , 'RESP2'   , 'RESP3'   , 
     1	'RESP4'   , 'RESP5'   , 'RESP6'   , 'RESP7'   , 'RESP8'   , 
     1	'RESP9'   , 'STLA'    , 'STLO'    , 'STEL'    , 'STDP'    , 
     1	'EVLA'    , 'EVLO'    , 'EVEL'    , 'EVDP'    , 'FHDR40'  , 
     1	'USER0'   , 'USER1'   , 'USER2'   , 'USER3'   , 'USER4'    /
	data (rstr(i),i=46,70)/
     1	'USER5'   ,'USER6'   , 'USER7'   , 'USER8'   , 'USER9'   , 
     1	'DIST'    ,'AZ'      , 'BAZ'     , 'GCARC'   , 'SB'      , 
     1	'SDELTA'  ,'DEPMEN'  , 'CMPAZ'   , 'CMPINC'  , 'XMINIMUM', 
     1	'XMAXIMUM','YMINIMUM', 'YMAXIMUM', 'ADJTM'   , 'FHDR65'  , 
     1	'FHDR66'  ,'FHDR67'  , 'FHDR68'  , 'FHDR69'  , 'FHDR70'   /
	data istr/
     1	'NZYEAR'  ,'NZJDAY'  ,'NZHOUR'  ,'NZMIN'   ,'NZSEC'   , 
     1	'NZMSEC'  ,'NVHDR'   ,'NINF'    ,'NHST'    ,'NPTS'    , 
     1	'NSNPTS'  ,'NSN'     ,'NXSIZE'  ,'NYSIZE'  ,'NHDR15'  , 
     1	'IFTYPE'  ,'IDEP'    ,'IZTYPE'  ,'IHDR4'   ,'IINST'   , 
     1	'ISTREG'  ,'IEVREG'  ,'IEVTYP'  ,'IQUAL'   ,'ISYNTH'  , 
     1	'IDHR11'  ,'IDHR12'  ,'IDHR13'  ,'IDHR14'  ,'IDHR15'  , 
     1	'IDHR16'  ,'IDHR17'  ,'IDHR18'  ,'IDHR19'  ,'IDHR20'  , 
     1	'LEVEN'   ,'LPSPOL'  ,'LOVROK'  ,'LCALDA'  ,'LHDR5'    /
     	data cstr/
     1	'KSTNM'   ,'KEVNM'   ,'KEVNMC'  ,'KHOLE'   , 
     1	'KO'      ,'KA'      ,'KT0'     ,'KT1'     , 
     1	'KT2'     ,'KT3'     ,'KT4'     ,'KT5'     , 
     1	'KT6'     ,'KT7'     ,'KT8'     ,'KT9'     , 
     1	'KF'      ,'KUSER0'  ,'KUSER1'  ,'KUSER2'  , 
     1	'KCMPNM'  ,'KNETWK'  ,'KDATRD'  ,'KINST'   
     1	/
c-----
c	formats
c-----
    9   format(' SAC Header with interval value and meaning')
   10	format(' REAL        INDEX  NAME',
     1		'        INT VALUE   REAL VALUE')
   11	format(13x,i5,2x,a8,e13.5,e13.5)
   20	format('INTEGER     INDEX  NAME',
     1		'        INT VALUE    INT VALUE')
   21	format(13x,i5,2x,a8,i13,i13)
   30	format(' ENNUMERATED INDEX  NAME',
     1		'        INT VALUE    ENU VALUE')
   31	format(13x,i5,2x,a8,i13,5x,a8)
   40	format(' LOGICAL     INDEX  NAME',
     1		'        INT VALUE    LOG VALUE')
   41	format(13x,i5,2x,a8,i13,5x,l8)
   50	format(' CHARACTER   INDEX  NAME',
     1		'        INT VALUE   CHAR VALUE')
   51	format(13x,i5,2x,a8,5x,a8,5x,a8)
c-----
c	output real header
c-----
	write(LOT,10)
	do 1000 i=1,70
		call getfhv(rstr(i),rval,nerr)
		write(LOT,11)i,rstr(i),rhdr(i), rval
 1000	continue
c-----
c	output integer header
c-----
	write(LOT,20)
C	do 2000 i=1,15
	do 2000 i=1,40
		call getnhv(istr(i),ival,nerr)
		write(LOT,21)i,istr(i),ihdr(i), ival
 2000	continue
c-----
c	output enumerated header
c-----
	write(LOT,30)
	do 3000 i=16,25
		call getihv(istr(i),cval,nerr)
		write(LOT,31)i,istr(i),ihdr(i), cval
 3000	continue
c-----
c	output logical header
c-----
	write(LOT,40)
	do 4000 i=36,40
		call getlhv(istr(i),lval,nerr)
		write(LOT,41)i,istr(i),ihdr(i), lval
 4000	continue
c-----
c	output character header
c-----
	write(LOT,50)
	do 5000 i=1,24
		call getkhv(cstr(i),cval,nett)
		write(LOT,51)i,cstr(i),chdr(i), cval
 5000	continue
	call getnhv('NPTS    ',npts,nerr)
	inc = npts/10
	inc = 1
	do 9000 i=1,npts,inc
		write(LOT,*)i,x(i)
 9000	continue
	write(LOT,*)npts,x(npts)
c-----
c	plot it
c-----
	npts = ihdr(10)
	return
	end

