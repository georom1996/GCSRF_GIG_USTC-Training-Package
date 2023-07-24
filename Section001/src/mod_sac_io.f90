MODULE sac_i_o

!------------------------------------------------------------------------------

! This module allows the reading and writing of sac files in either
! alphanumeric or binary format.  All the header info is included, 
! allowing for exact duplication of a read in/written out record. 

! 1) Compile with: 
!    f95 mod_sac_io.f90 "your_program" -kind=byte -o sac.x

! This typically produces the following warning (at least on my machine):
! Extension: mod_sac_io.f90, line 209: ALLOCATABLE dummy argument YARRAY
!           detected at ::@YARRAY
! Extension: mod_sac_io.f90, line 460: ALLOCATABLE dummy argument YARRAY
!           detected at ::@YARRAY
!
! This can be safely ignored.

! 2) Make sure to include the following line in your programs
!    USE sac_i_o

! Variables
! =======================================
! Sac header variables

! REALS: (32-bit or 4 bytes each)
! delta, depmin, depmax, scale, odelta
! b, e, o, a, internal1 
! t0, t1, t2, t3, t4 
! t5, t6, t7, t8, t9 
! f, resp0, resp1, resp2, resp3
! resp4, resp5, resp6, resp7, resp8
! resp9, stla, stlo, stel, stdp
! resp9, stla, stlo, stel, stdp
! evla, evlo, evel, evdp, mag
! user0, user1, user2, user3, user4
! user5, user6, user7, user8, user9
! dist, az, baz, gcarc, internal2
! internal3,  depmen, cmpaz, cmpinc, xminimum
! xmaximum,   yminimum,   ymaximum,   unused1, unused2
! unused3, unused4, unused5, unused6, unused7

! INTEGERS: (32-bit or 4 bytes each)
! nzyear, nzjday, nzhour, nzmin, nzsec 
! nzmsec, nvhdr, norid, nevid, npts
! internal4,   nwfid, nxsize, nysize, unused8
! iftype, idep, iztype, unused9, iinst
! istreg, ievreg, ievtyp, iqual, isynth
! imagtyp, imagsrc, unused10, unused11, unused12
! unused13, unused14, unused15, unused16, unused17

! LOGICALS: (32-bit or 4 bytes each)
! leven, lpspol, lovrok, lcalda, unused18      

! CHARACTERS: (64-bit or 8 bytes each)
! kstnm, kevnm*
! khole,ko,ka
! kt0,kt1,kt2
! kt3,kt4,kt5
! kt6,kt7,kt8
! kt9,kf,kuser0
! kuser1,kuser2,kcmpnm
! knetwk,kdatrd,kinst
! yarray
 
! *128 bit or 16 bytes

! Subroutines Included:
!=======================================
! rbsac: reads binary sacfile
! rasac: reads alphanumeric sacfile
! wbsac: writes binary sacfile
! wasac: writes alphanumeric sacfile


! To call the subroutines use (and replace rbsac w/ rasac, wasac, or wbsac):
!=======================================
! CALL rbsac(infile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,        &
! t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
! resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
! user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
! internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
! unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
! nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
! iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
! imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
! unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
! kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
! kinst,yarray)

! To undefine a header in your program:
!=======================================
! real_header_variable = rundef
! intg_header_variable = iundef
! logi_header_variable = ltrue -or- lfalse
! char_header_variable = kundef

!------------------------------------------------------------------------------

IMPLICIT NONE

! Header variables

! Select kind type (byte)
   INTEGER,PARAMETER :: k=4_4
! Real sac header variables
   REAL(k),PARAMETER :: rundef=-12345.0_4
   REAL(k) :: delta,depmin,depmax,scale,odelta
   REAL(k) :: b,e,o,a,internal1
   REAL(k) :: t0,t1,t2,t3,t4
   REAL(k) :: t5,t6,t7,t8,t9
   REAL(k) :: f,resp0,resp1,resp2,resp3
   REAL(k) :: resp4,resp5,resp6,resp7,resp8
   REAL(k) :: resp9,stla,stlo,stel,stdp
   REAL(k) :: evla,evlo,evel,evdp,mag
   REAL(k) :: user0,user1,user2,user3,user4
   REAL(k) :: user5,user6,user7,user8,user9
   REAL(k) :: dist,az,baz,gcarc,internal2
   REAL(k) :: internal3,depmen,cmpaz,cmpinc,xminimum
   REAL(k) :: xmaximum,yminimum,ymaximum,unused1,unused2
   REAL(k) :: unused3,unused4,unused5,unused6,unused7
! Integer sac header variables
   INTEGER(k),PARAMETER :: iundef=-12345_4
   INTEGER(k) :: nzyear,nzjday,nzhour,nzmin,nzsec
   INTEGER(k) :: nzmsec,nvhdr,norid,nevid,npts
   INTEGER(k) :: internal4,nwfid,nxsize,nysize,unused8
   INTEGER(k) :: iftype,idep,iztype,unused9,iinst
   INTEGER(k) :: istreg,ievreg,ievtyp,iqual,isynth
   INTEGER(k) :: imagtyp,imagsrc,unused10,unused11,unused12
   INTEGER(k) :: unused13,unused14,unused15,unused16,unused17
! Logical sac header variables
   INTEGER(k),PARAMETER :: ltrue=1_4,lfalse=0_4
   INTEGER(k) :: leven,lpspol,lovrok,lcalda,unused18
! Character sac header variables
   CHARACTER(LEN=8),PARAMETER :: kundef='-12345'
   CHARACTER(LEN=16) :: kevnm
   CHARACTER(LEN=8) :: kstnm
   CHARACTER(LEN=8) :: khole,ko,ka
   CHARACTER(LEN=8) :: kt0,kt1,kt2
   CHARACTER(LEN=8) :: kt3,kt4,kt5
   CHARACTER(LEN=8) :: kt6,kt7,kt8
   CHARACTER(LEN=8) :: kt9,kf,kuser0
   CHARACTER(LEN=8) :: kuser1,kuser2,kcmpnm
   CHARACTER(LEN=8) :: knetwk,kdatrd,kinst
! Data array
   REAL(k), ALLOCATABLE :: yarray(:)

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE rbsac(infile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,  &
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,yarray)

! Read in sac binary format file with ALL the header information
! ==============================
! Header values are as described on 
! the sac homepage http://www.llnl.gov/sac/

! Note: With NAGWare Fortran 95 Compiler
! must compile with -kind=byte; 
! has not been tested on non-Macintosh/Unix machines

   IMPLICIT NONE

   ! Real header variables
   INTEGER, PARAMETER :: k=4_4  !32 bit (4 byte) words
   REAL(k), INTENT(OUT):: delta,depmin,depmax,scale,odelta
   REAL(k), INTENT(OUT):: b,e,o,a,internal1
   REAL(k), INTENT(OUT):: t0,t1,t2,t3,t4
   REAL(k), INTENT(OUT):: t5,t6,t7,t8,t9
   REAL(k), INTENT(OUT):: f,resp0,resp1,resp2,resp3
   REAL(k), INTENT(OUT):: resp4,resp5,resp6,resp7,resp8
   REAL(k), INTENT(OUT):: resp9,stla,stlo,stel,stdp
   REAL(k), INTENT(OUT):: evla,evlo,evel,evdp,mag
   REAL(k), INTENT(OUT):: user0,user1,user2,user3,user4
   REAL(k), INTENT(OUT):: user5,user6,user7,user8,user9
   REAL(k), INTENT(OUT):: dist,az,baz,gcarc,internal2
   REAL(k), INTENT(OUT):: internal3,depmen,cmpaz,cmpinc,xminimum
   REAL(k), INTENT(OUT):: xmaximum,yminimum,ymaximum,unused1,unused2
   REAL(k), INTENT(OUT):: unused3,unused4,unused5,unused6,unused7

   ! Integer header variables
   INTEGER(k), INTENT(OUT):: nzyear,nzjday,nzhour,nzmin,nzsec
   INTEGER(k), INTENT(OUT):: nzmsec,nvhdr,norid,nevid,npts
   INTEGER(k), INTENT(OUT):: internal4,nwfid,nxsize,nysize,unused8
   INTEGER(k), INTENT(OUT):: iftype,idep,iztype,unused9,iinst
   INTEGER(k), INTENT(OUT):: istreg,ievreg,ievtyp,iqual,isynth
   INTEGER(k), INTENT(OUT):: imagtyp,imagsrc,unused10,unused11,unused12
   INTEGER(k), INTENT(OUT):: unused13,unused14,unused15,unused16,unused17

   ! Logical header variables (1=true 0=false)
   INTEGER(k), INTENT(OUT):: leven,lpspol,lovrok,lcalda,unused18

   ! Character header values
   CHARACTER(LEN=100), INTENT(IN) :: infile
   CHARACTER(LEN=16), INTENT(OUT):: kevnm
   CHARACTER(LEN=8), INTENT(OUT):: kstnm
   CHARACTER(LEN=8), INTENT(OUT):: khole,ko,ka
   CHARACTER(LEN=8), INTENT(OUT):: kt0,kt1,kt2
   CHARACTER(LEN=8), INTENT(OUT):: kt3,kt4,kt5
   CHARACTER(LEN=8), INTENT(OUT):: kt6,kt7,kt8
   CHARACTER(LEN=8), INTENT(OUT):: kt9,kf,kuser0
   CHARACTER(LEN=8), INTENT(OUT):: kuser1,kuser2,kcmpnm
   CHARACTER(LEN=8), INTENT(OUT):: knetwk,kdatrd,kinst
   ! Data array
   REAL(k), ALLOCATABLE,INTENT(OUT) :: yarray(:)

   ! Internal variables
   INTEGER(k) :: io,j,N,tot


      ! Read in float and integer data (each variable 4 units)
      OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=20)

      ! Position      Header Variables
      ! 1-20         DELTA DEPMIN DEPMAX SCALE ODELTA
      read(99,REC=1) delta, depmin, depmax, scale, odelta

      ! 21-40        B  E  O  A       INTERNAL
      read(99,REC=2) b, e, o, a, internal1 

      ! 41-60             T0      T1      T2      T3      T4
      read(99,REC=3) t0, t1, t2, t3, t4 

      ! 61-80             T5      T6      T7      T8      T9
      read(99,REC=4) t5, t6, t7, t8, t9 

      ! 81-100             F       RESP0   RESP1   RESP2   RESP3
      read(99,REC=5) f, resp0, resp1, resp2, resp3

      ! 101-120             RESP4   RESP5   RESP6   RESP7   RESP8 
      read(99,REC=6) resp4, resp5, resp6, resp7, resp8 

      ! 121-140             RESP9   STLA  STLO  STEL  STDP
      read(99,REC=7) resp9, stla, stlo, stel, stdp

      ! 141-160             EVLA  EVLO  EVEL    EVDP  MAG
      read(99,REC=8) evla, evlo, evel, evdp, mag

      ! 161-180             USER0   USER1   USER2   USER3   USER4
      read(99,REC=9) user0, user1, user2, user3, user4

      ! 181-200             USER5   USER6   USER7   USER8   USER9
      read(99,REC=10) user5, user6, user7, user8, user9

      ! 201-220             DIST    AZ      BAZ     GCARC   INTERNAL
      read(99,REC=11) dist, az, baz, gcarc, internal2

      ! 221-240             INTERNAL DEPMEN  CMPAZ  CMPINC  XMINIMUM
      read(99,REC=12) internal3,  depmen, cmpaz, cmpinc, xminimum

      ! 241-260             XMAXIMUM  YMINIMUM  YMAXIMUM  UNUSED  UNUSED
      read(99,REC=13) xmaximum,   yminimum,   ymaximum,   unused1, unused2

      ! 261-280             UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      read(99,REC=14) unused3, unused4, unused5, unused6, unused7
        
      !read integer header variables
      !----------------------------------------------------
      ! 281-300             NZYEAR  NZJDAY  NZHOUR  NZMIN  NZSEC
      read(99,REC=15) nzyear, nzjday, nzhour, nzmin, nzsec 

      ! 301-320             NZMSEC  NVHDR  NORID   NEVID   NPTS
      read(99,REC=16) nzmsec, nvhdr, norid, nevid, npts

      ! 321-340             INTERNAL  NWFID   NXSIZE  NYSIZE  UNUSED
      read(99,REC=17) internal4,   nwfid, nxsize, nysize, unused8

      ! 341-360             IFTYPE  IDEP    IZTYPE  UNUSED  IINST
      read(99,REC=18) iftype, idep, iztype, unused9, iinst

      ! 361-380             ISTREG  IEVREG  IEVTYP  IQUAL   ISYNTH
      read(99,REC=19) istreg, ievreg, ievtyp, iqual, isynth

      ! 381-400             IMAGTYP IMAGSRC UNUSED  UNUSED  UNUSED
      read(99,REC=20) imagtyp, imagsrc, unused10, unused11, unused12

      ! 401-420             UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      read(99,REC=21) unused13, unused14, unused15, unused16, unused17

      !read logical header variables
      !----------------------------------------------------
      ! 421-440             LEVEN  LPSPOL  LOVROK LCALDA UNUSED
      read(99,REC=22) leven, lpspol, lovrok, lcalda, unused18
      CLOSE(99)

      ! Character variable (size=8)*   *kevnm=16
      ! Here it gets tricky; we have to read in each k-header separately
      ! due to their different size... so I deviate from a straight up
      ! copy of the alpha numeric format
      
      !read character header variables
      !----------------------------------------------------
      ! 440-448             KSTNM
      OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=8)
      read(99,REC=56) kstnm
      CLOSE(99)

      ! 449-464             KEVNM
      OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=16)
      read(99,REC=29) kevnm

      ! 465-480             KHOLE   KO    
      read(99,REC=30) khole, ko
      CLOSE(99)

      ! 480-504             KA   KT0     KT1     
      OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=24)
      read(99,REC=21) ka, kt0, kt1

      ! 504-528             KT2     KT3     KT4
      read(99,REC=22) kt2, kt3, kt4
      ! 529-552              KT5     KT6     KT7
      read(99,REC=23) kt5, kt6, kt7

      ! 553-576             KT8  KT9   KF     
      read(99,REC=24) kt8, kt9, kf 

      ! 577-600             KUSER0 KUSER1  KUSER2  
      read(99,REC=25) kuser0, kuser1, kuser2

      ! 601-624             KCMPNM KNETWK  KDATRD  
      read(99,REC=26) kcmpnm, knetwk, kdatrd
      CLOSE(99)

      ! 625-632             KINST
      OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=8)
      read(99,REC=79) kinst
      CLOSE(99)

      !read data
      ! 633 ->  npts*4
      !----------------------------------------------------
      ALLOCATE (yarray(npts))


      ! Efficiency counts; I've worked out a way of reading in multiple
      ! binary values at once

      ! 1) If there are less than 158 points, then a simple do loop is 
      ! best
      ! 2) If there are more than 158 points, you want to read in blocks
      ! of 158 pts each (Direct access requires you keep track of where
      ! you are in the file; multiples of 158 are good as there are 158
      ! header variables, so you count them as record 1).
      ! However, you also need to read in the rest of the points as it
      ! is unlikely files will always have a multiple of 158 for npts

      ! Case 1: less than 158 pts in the sacfile
      IF (npts < 158_4) THEN
!WRITE(*,*) "npts < 158"

        OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=4)
        DO j=159_4,(159_4+npts)

           READ(99,REC=j,iostat=io) yarray(j-158_4)
           IF (io /= 0) EXIT

        END DO
        CLOSE(99)
         
      ! Case 2: 
      ELSE IF (npts >= 158_4) THEN
!WRITE(*,*) "npts > 158"
       
        OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=632)

        N=int(real(npts)/158_4)+1_4

!WRITE(*,*) N,npts

        ! Start reading in the large record blocks (158*4 blocks each)
        DO j=2_4,N

           READ(99,REC=j,iostat=io) yarray(158_4*(j-2_4)+1_4:158_4*(j-1_4))
!WRITE(*,*) io
           IF (io /= 0) EXIT
        END DO
        CLOSE(99)

        tot=npts-158_4*int(real(npts)/158_4)

!WRITE(*,*) tot

        ! Case 3: Read in the rest
        OPEN(UNIT=99,FILE=infile,FORM='unformatted',ACCESS='DIRECT',RECL=4)
        DO j=1,tot
           READ(99,REC=((158_4*N)+j),iostat=io) yarray(158_4*(N-1_4)+j)
! WRITE(*,*) 158_4*(N)+j,yarray(158_4*(N)+j),io
!           IF (io /= 0) EXIT
        END DO
        CLOSE(99)

      END IF


! For Debugging:
!      write(*,201) delta, depmin, depmax, scale, odelta
!      write(*,201) b, e, o, a, internal1 
!      write(*,201) t0, t1, t2, t3, t4 
!      write(*,201) t5, t6, t7, t8, t9 
!      write(*,201) f, resp0, resp1, resp2, resp3
!      write(*,201) resp4, resp5, resp6, resp7, resp8
!      write(*,201) resp9, stla, stlo, stel, stdp
!      write(*,201) resp9, stla, stlo, stel, stdp
!      write(*,201) evla, evlo, evel, evdp, mag
!      write(*,201) user0, user1, user2, user3, user4
!      write(*,201) user5, user6, user7, user8, user9
!      write(*,201) dist, az, baz, gcarc, internal2
!      write(*,201) internal3,  depmen, cmpaz, cmpinc, xminimum
!      write(*,201) xmaximum,   yminimum,   ymaximum,   unused1, unused2
!      write(*,201) unused3, unused4, unused5, unused6, unused7
!      write(*,202) nzyear, nzjday, nzhour, nzmin, nzsec 
!      write(*,202) nzmsec, nvhdr, norid, nevid, npts
!      write(*,202) internal4,   nwfid, nxsize, nysize, unused8
!      write(*,202) iftype, idep, iztype, unused9, iinst
!      write(*,202) istreg, ievreg, ievtyp, iqual, isynth
!      write(*,202) imagtyp, imagsrc, unused10, unused11, unused12
!      write(*,202) unused13, unused14, unused15, unused16, unused17
!      write(*,202) leven, lpspol, lovrok, lcalda, unused18      
!      write(*,203) kstnm, kevnm
!      write(*,204) khole,ko,ka
!      write(*,204) kt0,kt1,kt2
!      write(*,204) kt3,kt4,kt5
!      write(*,204) kt6,kt7,kt8
!      write(*,204) kt9,kf,kuser0
!      write(*,204) kuser1,kuser2,kcmpnm
!      write(*,204) knetwk,kdatrd,kinst
!      write(*,205) yarray

!    201 FORMAT(G15.7,G15.7,G15.7,G15.7,G15.7)
!    202 FORMAT(I10,I10,I10,I10,I10)
!    203 FORMAT(A8,A16)
!    204 FORMAT(A8,A8,A8)
!    205 FORMAT(G15.7,G15.7,G15.7,G15.7,G15.7)

  END SUBROUTINE rbsac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE rasac(infile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,&
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,yarray)

! Read in sac alphanumeric format file with ALL the header information
! ==============================
! Header values are as described on 
! the sac homepage http://www.llnl.gov/sac/

! Note: With NAGWare Fortran 95 Compiler
! must compile with -kind=byte; 
! has not been tested on non-Macintosh/Unix machines

   IMPLICIT NONE

   ! Real header variables
   REAL, INTENT(OUT):: delta,depmin,depmax,scale,odelta
   REAL, INTENT(OUT):: b,e,o,a,internal1
   REAL, INTENT(OUT):: t0,t1,t2,t3,t4
   REAL, INTENT(OUT):: t5,t6,t7,t8,t9
   REAL, INTENT(OUT):: f,resp0,resp1,resp2,resp3
   REAL, INTENT(OUT):: resp4,resp5,resp6,resp7,resp8
   REAL, INTENT(OUT):: resp9,stla,stlo,stel,stdp
   REAL, INTENT(OUT):: evla,evlo,evel,evdp,mag
   REAL, INTENT(OUT):: user0,user1,user2,user3,user4
   REAL, INTENT(OUT):: user5,user6,user7,user8,user9
   REAL, INTENT(OUT):: dist,az,baz,gcarc,internal2
   REAL, INTENT(OUT):: internal3,depmen,cmpaz,cmpinc,xminimum
   REAL, INTENT(OUT):: xmaximum,yminimum,ymaximum,unused1,unused2
   REAL, INTENT(OUT):: unused3,unused4,unused5,unused6,unused7

   ! integer header variables
   INTEGER, INTENT(OUT):: nzyear,nzjday,nzhour,nzmin,nzsec
   INTEGER, INTENT(OUT):: nzmsec,nvhdr,norid,nevid,npts
   INTEGER, INTENT(OUT):: internal4,nwfid,nxsize,nysize,unused8
   INTEGER, INTENT(OUT):: iftype,idep,iztype,unused9,iinst
   INTEGER, INTENT(OUT):: istreg,ievreg,ievtyp,iqual,isynth
   INTEGER, INTENT(OUT):: imagtyp,imagsrc,unused10,unused11,unused12
   INTEGER, INTENT(OUT):: unused13,unused14,unused15,unused16,unused17

   ! logical header variables (1=true 0=false)
   INTEGER, INTENT(OUT):: leven,lpspol,lovrok,lcalda,unused18

   ! character header values
   CHARACTER(LEN=90), INTENT(IN) :: infile
   CHARACTER(LEN=16), INTENT(OUT):: kevnm
   CHARACTER(LEN=8), INTENT(OUT):: kstnm
   CHARACTER(LEN=8), INTENT(OUT):: khole,ko,ka
   CHARACTER(LEN=8), INTENT(OUT):: kt0,kt1,kt2
   CHARACTER(LEN=8), INTENT(OUT):: kt3,kt4,kt5
   CHARACTER(LEN=8), INTENT(OUT):: kt6,kt7,kt8
   CHARACTER(LEN=8), INTENT(OUT):: kt9,kf,kuser0
   CHARACTER(LEN=8), INTENT(OUT):: kuser1,kuser2,kcmpnm
   CHARACTER(LEN=8), INTENT(OUT):: knetwk,kdatrd,kinst

   ! Data array
   REAL, ALLOCATABLE,INTENT(OUT) :: yarray(:)

   ! status='replace' removed
   OPEN(UNIT=777,FILE=infile)
    
      !read in real header variables
      !----------------------------------------------------
      !              DELTA  DEPMIN  DEPMAX  SCALE   ODELTA
      read(777,201) delta, depmin, depmax, scale, odelta

      !              B  E  O  A       INTERNAL
      read(777,201) b, e, o, a, internal1 

      !              T0      T1      T2      T3      T4
      read(777,201) t0, t1, t2, t3, t4 

      !              T5      T6      T7      T8      T9
      read(777,201) t5, t6, t7, t8, t9 

      !              F       RESP0   RESP1   RESP2   RESP3
      read(777,201) f, resp0, resp1, resp2, resp3

      !              RESP4   RESP5   RESP6   RESP7   RESP8 
      read(777,201) resp4, resp5, resp6, resp7, resp8 

      !              RESP9   STLA  STLO  STEL  STDP
      read(777,201) resp9, stla, stlo, stel, stdp

      !              EVLA  EVLO  EVEL    EVDP  MAG
      read(777,201) evla, evlo, evel, evdp, mag

      !              USER0   USER1   USER2   USER3   USER4
      read(777,201) user0, user1, user2, user3, user4

      !              USER5   USER6   USER7   USER8   USER9
      read(777,201) user5, user6, user7, user8, user9

      !              DIST    AZ      BAZ     GCARC   INTERNAL
      read(777,201) dist, az, baz, gcarc, internal2

      !              INTERNAL DEPMEN  CMPAZ  CMPINC  XMINIMUM
      read(777,201) internal3,  depmen, cmpaz, cmpinc, xminimum

      !              XMAXIMUM  YMINIMUM  YMAXIMUM  UNUSED  UNUSED
      read(777,201) xmaximum,   yminimum,   ymaximum,   unused1, unused2

      !              UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      read(777,201) unused3, unused4, unused5, unused6, unused7

        
      !read integer header variables
      !----------------------------------------------------
      !              NZYEAR  NZJDAY  NZHOUR  NZMIN  NZSEC
      read(777,202) nzyear, nzjday, nzhour, nzmin, nzsec 

      !              NZMSEC  NVHDR  NORID   NEVID   NPTS
      read(777,202) nzmsec, nvhdr, norid, nevid, npts

      !              INTERNAL  NWFID   NXSIZE  NYSIZE  UNUSED
      read(777,202) internal4,   nwfid, nxsize, nysize, unused8

      !              IFTYPE  IDEP    IZTYPE  UNUSED  IINST
      read(777,202) iftype, idep, iztype, unused9, iinst

      !              ISTREG  IEVREG  IEVTYP  IQUAL   ISYNTH
      read(777,202) istreg, ievreg, ievtyp, iqual, isynth

      !              IMAGTYP IMAGSRC UNUSED  UNUSED  UNUSED
      read(777,202) imagtyp, imagsrc, unused10, unused11, unused12

      !              UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      read(777,202) unused13, unused14, unused15, unused16, unused17

      !read logical header variables
      !----------------------------------------------------
      !              LEVEN  LPSPOL  LOVROK LCALDA UNUSED
      read(777,202) leven, lpspol, lovrok, lcalda, unused18

      !read character header variables
      !----------------------------------------------------
      !              KSTNM  KEVNM
      read(777,203) kstnm, kevnm

      !              KHOLE   KO      KA
      read(777,204) khole, ko, ka

      !              KT0     KT1     KT2 
      read(777,204) kt0, kt1, kt2

      !              KT3     KT4     KT5
      read(777,204) kt3, kt4, kt5

      !              KT6     KT7     KT8
      read(777,204) kt6, kt7, kt8

      !              KT9     KF      KUSER0
      read(777,204) kt9, kf, kuser0

      !              KUSER1  KUSER2  KCMPNM
      read(777,204) kuser1, kuser2, kcmpnm

      !              KNETWK  KDATRD  KINST
      read(777,204) knetwk, kdatrd, kinst


      !read data
      !----------------------------------------------------
      ALLOCATE (yarray(npts))

      read(777,205) yarray

    CLOSE(777)

    201 FORMAT(G15.7,G15.7,G15.7,G15.7,G15.7)
    202 FORMAT(I10,I10,I10,I10,I10)
    203 FORMAT(A8,A16)
    204 FORMAT(A8,A8,A8)
    205 FORMAT(G15.7,G15.7,G15.7,G15.7,G15.7)

  END SUBROUTINE rasac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE wbsac(outfile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,&
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,yarray)

! Write out sac binary format file with ALL the header information
! ==============================
! Header values are as described on 
! the sac homepage http://www.llnl.gov/sac/

! Note: With NAGWare Fortran 95 Compiler
! must compile with -kind=byte; 
! has not been tested on non-Macintosh/Unix machines

   IMPLICIT NONE

   ! Real header variables
   INTEGER, PARAMETER :: k=4_4  !32 bit (4 byte) words
   REAL(k), INTENT(IN):: delta,depmin,depmax,scale,odelta
   REAL(k), INTENT(IN):: b,e,o,a,internal1
   REAL(k), INTENT(IN):: t0,t1,t2,t3,t4
   REAL(k), INTENT(IN):: t5,t6,t7,t8,t9
   REAL(k), INTENT(IN):: f,resp0,resp1,resp2,resp3
   REAL(k), INTENT(IN):: resp4,resp5,resp6,resp7,resp8
   REAL(k), INTENT(IN):: resp9,stla,stlo,stel,stdp
   REAL(k), INTENT(IN):: evla,evlo,evel,evdp,mag
   REAL(k), INTENT(IN):: user0,user1,user2,user3,user4
   REAL(k), INTENT(IN):: user5,user6,user7,user8,user9
   REAL(k), INTENT(IN):: dist,az,baz,gcarc,internal2
   REAL(k), INTENT(IN):: internal3,depmen,cmpaz,cmpinc,xminimum
   REAL(k), INTENT(IN):: xmaximum,yminimum,ymaximum,unused1,unused2
   REAL(k), INTENT(IN):: unused3,unused4,unused5,unused6,unused7

   ! integer header variables
   INTEGER(k), INTENT(IN):: nzyear,nzjday,nzhour,nzmin,nzsec
   INTEGER(k), INTENT(IN):: nzmsec,nvhdr,norid,nevid,npts
   INTEGER(k), INTENT(IN):: internal4,nwfid,nxsize,nysize,unused8
   INTEGER(k), INTENT(IN):: iftype,idep,iztype,unused9,iinst
   INTEGER(k), INTENT(IN):: istreg,ievreg,ievtyp,iqual,isynth
   INTEGER(k), INTENT(IN):: imagtyp,imagsrc,unused10,unused11,unused12
   INTEGER(k), INTENT(IN):: unused13,unused14,unused15,unused16,unused17

   ! logical header variables (1=true 0=false)
   INTEGER(k), INTENT(IN):: leven,lpspol,lovrok,lcalda,unused18

   ! character header values
   CHARACTER(LEN=100), INTENT(IN) :: outfile
   CHARACTER(LEN=16), INTENT(IN):: kevnm
   CHARACTER(LEN=8), INTENT(IN):: kstnm
   CHARACTER(LEN=8), INTENT(IN):: khole,ko,ka
   CHARACTER(LEN=8), INTENT(IN):: kt0,kt1,kt2
   CHARACTER(LEN=8), INTENT(IN):: kt3,kt4,kt5
   CHARACTER(LEN=8), INTENT(IN):: kt6,kt7,kt8
   CHARACTER(LEN=8), INTENT(IN):: kt9,kf,kuser0
   CHARACTER(LEN=8), INTENT(IN):: kuser1,kuser2,kcmpnm
   CHARACTER(LEN=8), INTENT(IN):: knetwk,kdatrd,kinst
   ! Data array
   REAL(k),INTENT(IN) :: yarray(npts)

   ! Internal variables
   INTEGER(k) :: io,j,N,tot

      ! Write out float and integer data (4 units/variable)
      OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=20)

      ! Position      Header Variables
      ! 1-20         DELTA DEPMIN DEPMAX SCALE ODELTA
      write(20,REC=1) delta, depmin, depmax, scale, odelta

      ! 21-40        B  E  O  A       INTERNAL
      write(20,REC=2) b, e, o, a, internal1 

      ! 41-60             T0      T1      T2      T3      T4
      write(20,REC=3) t0, t1, t2, t3, t4 

      ! 61-80             T5      T6      T7      T8      T9
      write(20,REC=4) t5, t6, t7, t8, t9 

      ! 81-100             F       RESP0   RESP1   RESP2   RESP3
      write(20,REC=5) f, resp0, resp1, resp2, resp3

      ! 101-120             RESP4   RESP5   RESP6   RESP7   RESP8 
      write(20,REC=6) resp4, resp5, resp6, resp7, resp8 

      ! 121-140             RESP9   STLA  STLO  STEL  STDP
      write(20,REC=7) resp9, stla, stlo, stel, stdp

      ! 141-160             EVLA  EVLO  EVEL    EVDP  MAG
      write(20,REC=8) evla, evlo, evel, evdp, mag

      ! 161-180             USER0   USER1   USER2   USER3   USER4
      write(20,REC=9) user0, user1, user2, user3, user4

      ! 181-200             USER5   USER6   USER7   USER8   USER9
      write(20,REC=10) user5, user6, user7, user8, user9

      ! 201-220             DIST    AZ      BAZ     GCARC   INTERNAL
      write(20,REC=11) dist, az, baz, gcarc, internal2

      ! 221-240             INTERNAL DEPMEN  CMPAZ  CMPINC  XMINIMUM
      write(20,REC=12) internal3,  depmen, cmpaz, cmpinc, xminimum

      ! 241-260             XMAXIMUM  YMINIMUM  YMAXIMUM  UNUSED  UNUSED
      write(20,REC=13) xmaximum,   yminimum,   ymaximum,   unused1, unused2

      ! 261-280             UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      write(20,REC=14) unused3, unused4, unused5, unused6, unused7
        
      !write integer header variables
      !----------------------------------------------------
      ! 281-300             NZYEAR  NZJDAY  NZHOUR  NZMIN  NZSEC
      write(20,REC=15) nzyear, nzjday, nzhour, nzmin, nzsec 

      ! 301-320             NZMSEC  NVHDR  NORID   NEVID   NPTS
      write(20,REC=16) nzmsec, nvhdr, norid, nevid, npts

      ! 321-340             INTERNAL  NWFID   NXSIZE  NYSIZE  UNUSED
      write(20,REC=17) internal4,   nwfid, nxsize, nysize, unused8

      ! 341-360             IFTYPE  IDEP    IZTYPE  UNUSED  IINST
      write(20,REC=18) iftype, idep, iztype, unused9, iinst

      ! 361-380             ISTREG  IEVREG  IEVTYP  IQUAL   ISYNTH
      write(20,REC=19) istreg, ievreg, ievtyp, iqual, isynth

      ! 381-400             IMAGTYP IMAGSRC UNUSED  UNUSED  UNUSED
      write(20,REC=20) imagtyp, imagsrc, unused10, unused11, unused12

      ! 401-420             UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      write(20,REC=21) unused13, unused14, unused15, unused16, unused17

      !write logical header variables
      !----------------------------------------------------
      ! 421-440             LEVEN  LPSPOL  LOVROK LCALDA UNUSED
      write(20,REC=22) leven, lpspol, lovrok, lcalda, unused18
      CLOSE(20)

      ! Character variable (size=8)*   *kevnm=16
      ! Here it gets tricky; we have to write out each k-header separately
      ! due to their different size... so I deviate from a straight up
      ! copy of the alpha numeric format, and do it my own way
      
      !write character header variables
      !----------------------------------------------------
      ! 440-448             KSTNM
      OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=8)
      write(20,REC=56) kstnm
      CLOSE(20)
      ! 449-464             KEVNM
      OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=16)
      write(20,REC=29) kevnm

      ! 465-480             KHOLE   KO    
      write(20,REC=30) khole, ko
      CLOSE(20)

      ! 480-504             KA   KT0     KT1     
      OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=24)
      write(20,REC=21) ka, kt0, kt1

      ! 504-528             KT2     KT3     KT4
      write(20,REC=22) kt2, kt3, kt4

      ! 529-552              KT5     KT6     KT7
      write(20,REC=23) kt5, kt6, kt7

      ! 553-576             KT8  KT9   KF     
      write(20,REC=24) kt8, kt9, kf 

      ! 577-600             KUSER0 KUSER1  KUSER2  
      write(20,REC=25) kuser0, kuser1, kuser2

      ! 601-624             KCMPNM KNETWK  KDATRD  
      write(20,REC=26) kcmpnm, knetwk, kdatrd
      CLOSE(20)

      ! 625-632             KINST
      OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=8)
      write(20,REC=79) kinst
      CLOSE(20)

      ! write data
      ! 633 ->  npts*4
      !----------------------------------------------------

!      OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=4)
!
!      DO j=159_4,(159_4+npts)
!
!         WRITE(20,REC=j,iostat=io) yarray(j-158_4)
!         IF (io /= 0) EXIT
!
!      END DO
!
!      CLOSE(20)

      ! Efficiency counts; I've worked out a way of reading in multiple
      ! binary values at once

      ! 1) If there are less than 158 points, then a simple do loop is 
      ! best
      ! 2) If there are more than 158 points, you want to read in blocks
      ! of 158 pts each (Direct access requires you keep track of where
      ! you are in the file; multiples of 158 are good as there are 158
      ! header variables, so you count them as record 1).
      ! However, you also need to read in the rest of the points as it
      ! is unlikely files will always have a multiple of 158 for npts

      ! Case 1: less than 158 pts in the sacfile
      IF (npts < 158_4) THEN

        OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=4)
        DO j=159_4,(159_4+npts)

           WRITE(20,REC=j,iostat=io) yarray(j-158_4)
           IF (io /= 0) EXIT

        END DO
        CLOSE(20)
         
      ! Case 2: 
      ELSE IF (npts >= 158_4) THEN
        
        OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=632)

        N=int(real(npts)/158_4)+1_4

        ! Start reading in the large record blocks (158*4 blocks each)
        DO j=2_4,N

           WRITE(20,REC=j,iostat=io) yarray(158_4*(j-2_4)+1_4:158_4*(j-1_4))
           IF (io /= 0) EXIT
        END DO
        CLOSE(20)

        tot=npts-158_4*int(real(npts)/158_4)

        ! Case 3: Read in the rest
        OPEN(UNIT=20,FILE=outfile,FORM='unformatted',ACCESS='DIRECT',RECL=4)
        DO j=1,tot
           WRITE(20,REC=((158_4*N)+j),iostat=io) yarray(158_4*(N-1_4)+j)
!           WRITE(20,REC=((158_4*N)+j),iostat=io) yarray(158_4*(N-1)+j:npts)
           IF (io /= 0) EXIT
        END DO

      END IF

  END SUBROUTINE wbsac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE wasac(outfile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,&
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,    &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,    &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2,  &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,    &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,   &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,    &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,   &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,       &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1, &
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd, &
kinst,yarray)

! Write out sac alphanumeric format file with ALL the header information
! ==============================
! Header values are as described on 
! the sac homepage http://www.llnl.gov/sac/

! Note: With NAGWare Fortran 95 Compiler
! must compile with -kind=byte; 
! has not been tested on non-Macintosh/Unix machines

   IMPLICIT NONE

   ! Real header variables
   REAL, INTENT(IN):: delta,depmin,depmax,scale,odelta
   REAL, INTENT(IN):: b,e,o,a,internal1
   REAL, INTENT(IN):: t0,t1,t2,t3,t4
   REAL, INTENT(IN):: t5,t6,t7,t8,t9
   REAL, INTENT(IN):: f,resp0,resp1,resp2,resp3
   REAL, INTENT(IN):: resp4,resp5,resp6,resp7,resp8
   REAL, INTENT(IN):: resp9,stla,stlo,stel,stdp
   REAL, INTENT(IN):: evla,evlo,evel,evdp,mag
   REAL, INTENT(IN):: user0,user1,user2,user3,user4
   REAL, INTENT(IN):: user5,user6,user7,user8,user9
   REAL, INTENT(IN):: dist,az,baz,gcarc,internal2
   REAL, INTENT(IN):: internal3,depmen,cmpaz,cmpinc,xminimum
   REAL, INTENT(IN):: xmaximum,yminimum,ymaximum,unused1,unused2
   REAL, INTENT(IN):: unused3,unused4,unused5,unused6,unused7

   ! integer header variables
   INTEGER, INTENT(IN):: nzyear,nzjday,nzhour,nzmin,nzsec
   INTEGER, INTENT(IN):: nzmsec,nvhdr,norid,nevid,npts
   INTEGER, INTENT(IN):: internal4,nwfid,nxsize,nysize,unused8
   INTEGER, INTENT(IN):: iftype,idep,iztype,unused9,iinst
   INTEGER, INTENT(IN):: istreg,ievreg,ievtyp,iqual,isynth
   INTEGER, INTENT(IN):: imagtyp,imagsrc,unused10,unused11,unused12
   INTEGER, INTENT(IN):: unused13,unused14,unused15,unused16,unused17

   ! logical header variables (1=true 0=false)
   INTEGER, INTENT(IN):: leven,lpspol,lovrok,lcalda,unused18

   ! character header values
   CHARACTER(LEN=90), INTENT(IN) :: outfile
   CHARACTER(LEN=16), INTENT(IN):: kevnm
   CHARACTER(LEN=8), INTENT(IN):: kstnm
   CHARACTER(LEN=8), INTENT(IN):: khole,ko,ka
   CHARACTER(LEN=8), INTENT(IN):: kt0,kt1,kt2
   CHARACTER(LEN=8), INTENT(IN):: kt3,kt4,kt5
   CHARACTER(LEN=8), INTENT(IN):: kt6,kt7,kt8
   CHARACTER(LEN=8), INTENT(IN):: kt9,kf,kuser0
   CHARACTER(LEN=8), INTENT(IN):: kuser1,kuser2,kcmpnm
   CHARACTER(LEN=8), INTENT(IN):: knetwk,kdatrd,kinst
   
   ! Data array
   REAL, INTENT(IN) :: yarray(npts)
    
   ! Open the new sac file
   OPEN(UNIT=777,FILE=outfile,status='replace')
    
      !read in real header variables
      !----------------------------------------------------
      !              DELTA  DEPMIN  DEPMAX  SCALE   ODELTA
      write(777,201) delta, depmin, depmax, scale, odelta

      !              B  E  O  A       INTERNAL
      write(777,201) b, e, o, a, internal1 

      !              T0      T1      T2      T3      T4
      write(777,201) t0, t1, t2, t3, t4 

      !              T5      T6      T7      T8      T9
      write(777,201) t5, t6, t7, t8, t9 

      !              F       RESP0   RESP1   RESP2   RESP3
      write(777,201) f, resp0, resp1, resp2, resp3

      !              RESP4   RESP5   RESP6   RESP7   RESP8 
      write(777,201) resp4, resp5, resp6, resp7, resp8 

      !              RESP9   STLA  STLO  STEL  STDP
      write(777,201) resp9, stla, stlo, stel, stdp

      !              EVLA  EVLO  EVEL    EVDP  MAG
      write(777,201) evla, evlo, evel, evdp, mag

      !              USER0   USER1   USER2   USER3   USER4
      write(777,201) user0, user1, user2, user3, user4

      !              USER5   USER6   USER7   USER8   USER9
      write(777,201) user5, user6, user7, user8, user9

      !              DIST    AZ      BAZ     GCARC   INTERNAL
      write(777,201) dist, az, baz, gcarc, internal2

      !              INTERNAL DEPMEN  CMPAZ  CMPINC  XMINIMUM
      write(777,201) internal3,  depmen, cmpaz, cmpinc, xminimum

      !              XMAXIMUM  YMINIMUM  YMAXIMUM  UNUSED  UNUSED
      write(777,201) xmaximum,   yminimum,   ymaximum,   unused1, unused2

      !              UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      write(777,201) unused3, unused4, unused5, unused6, unused7

        
      !write integer header variables
      !----------------------------------------------------
      !              NZYEAR  NZJDAY  NZHOUR  NZMIN  NZSEC
      write(777,202) nzyear, nzjday, nzhour, nzmin, nzsec 

      !              NZMSEC  NVHDR  NORID   NEVID   NPTS
      write(777,202) nzmsec, nvhdr, norid, nevid, npts

      !              INTERNAL  NWFID   NXSIZE  NYSIZE  UNUSED
      write(777,202) internal4,   nwfid, nxsize, nysize, unused8

      !              IFTYPE  IDEP    IZTYPE  UNUSED  IINST
      write(777,202) iftype, idep, iztype, unused9, iinst

      !              ISTREG  IEVREG  IEVTYP  IQUAL   ISYNTH
      write(777,202) istreg, ievreg, ievtyp, iqual, isynth

      !              IMAGTYP IMAGSRC UNUSED  UNUSED  UNUSED
      write(777,202) imagtyp, imagsrc, unused10, unused11, unused12

      !              UNUSED  UNUSED  UNUSED  UNUSED  UNUSED
      write(777,202) unused13, unused14, unused15, unused16, unused17

      !write logical header variables
      !----------------------------------------------------
      !              LEVEN  LPSPOL  LOVROK LCALDA UNUSED
      write(777,202) leven, lpspol, lovrok, lcalda, unused18

      !write character header variables
      !----------------------------------------------------
      !              KSTNM  KEVNM
      write(777,203) kstnm, kevnm

      !              KHOLE   KO      KA
      write(777,204) khole, ko, ka

      !              KT0     KT1     KT2 
      write(777,204) kt0, kt1, kt2

      !              KT3     KT4     KT5
      write(777,204) kt3, kt4, kt5

      !              KT6     KT7     KT8
      write(777,204) kt6, kt7, kt8

      !              KT9     KF      KUSER0
      write(777,204) kt9, kf, kuser0

      !              KUSER1  KUSER2  KCMPNM
      write(777,204) kuser1, kuser2, kcmpnm

      !              KNETWK  KDATRD  KINST
      write(777,204) knetwk, kdatrd, kinst


      !write data
      !----------------------------------------------------

      write(777,205) yarray

    CLOSE(777)

    201 FORMAT(G15.7,G15.7,G15.7,G15.7,G15.7)
    202 FORMAT(I10,I10,I10,I10,I10)
    203 FORMAT(A8,A16)
    204 FORMAT(A8,A8,A8)
    205 FORMAT(G15.7,G15.7,G15.7,G15.7,G15.7)

  END SUBROUTINE wasac

  SUBROUTINE initsac(outfile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,  &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,    &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,    &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2,  &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,    &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,   &
    nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,    &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,   &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,       &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1, &
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd, &
    kinst,yarray)

   ! Real header variables
   INTEGER, PARAMETER :: k=4_4  !32 bit (4 byte) words
   REAL(k), INTENT(OUT):: delta,depmin,depmax,scale,odelta
   REAL(k), INTENT(OUT):: b,e,o,a,internal1
   REAL(k), INTENT(OUT):: t0,t1,t2,t3,t4
   REAL(k), INTENT(OUT):: t5,t6,t7,t8,t9
   REAL(k), INTENT(OUT):: f,resp0,resp1,resp2,resp3
   REAL(k), INTENT(OUT):: resp4,resp5,resp6,resp7,resp8
   REAL(k), INTENT(OUT):: resp9,stla,stlo,stel,stdp
   REAL(k), INTENT(OUT):: evla,evlo,evel,evdp,mag
   REAL(k), INTENT(OUT):: user0,user1,user2,user3,user4
   REAL(k), INTENT(OUT):: user5,user6,user7,user8,user9
   REAL(k), INTENT(OUT):: dist,az,baz,gcarc,internal2
   REAL(k), INTENT(OUT):: internal3,depmen,cmpaz,cmpinc,xminimum
   REAL(k), INTENT(OUT):: xmaximum,yminimum,ymaximum,unused1,unused2
   REAL(k), INTENT(OUT):: unused3,unused4,unused5,unused6,unused7
   REAL(k), PARAMETER  :: rundef=-12345.0

   ! Integer header variables
   INTEGER(k), INTENT(OUT):: nzyear,nzjday,nzhour,nzmin,nzsec
   INTEGER(k), INTENT(OUT):: nzmsec,nvhdr,norid,nevid,npts
   INTEGER(k), INTENT(OUT):: internal4,nwfid,nxsize,nysize,unused8
   INTEGER(k), INTENT(OUT):: iftype,idep,iztype,unused9,iinst
   INTEGER(k), INTENT(OUT):: istreg,ievreg,ievtyp,iqual,isynth
   INTEGER(k), INTENT(OUT):: imagtyp,imagsrc,unused10,unused11,unused12
   INTEGER(k), INTENT(OUT):: unused13,unused14,unused15,unused16,unused17
   INTEGER(k), PARAMETER  :: iundef=-12345

   ! Logical header variables (1=true 0=false)
   INTEGER(k), INTENT(OUT):: leven,lpspol,lovrok,lcalda,unused18

   ! Character header values
   CHARACTER(LEN=100), INTENT(OUT) :: outfile
   CHARACTER(LEN=16), INTENT(OUT):: kevnm
   CHARACTER(LEN=8), INTENT(OUT):: kstnm
   CHARACTER(LEN=8), INTENT(OUT):: khole,ko,ka
   CHARACTER(LEN=8), INTENT(OUT):: kt0,kt1,kt2
   CHARACTER(LEN=8), INTENT(OUT):: kt3,kt4,kt5
   CHARACTER(LEN=8), INTENT(OUT):: kt6,kt7,kt8
   CHARACTER(LEN=8), INTENT(OUT):: kt9,kf,kuser0
   CHARACTER(LEN=8), INTENT(OUT):: kuser1,kuser2,kcmpnm
   CHARACTER(LEN=8), INTENT(OUT):: knetwk,kdatrd,kinst
   CHARACTER(LEN=8)             :: kundef='-12345'
   ! Data array
   REAL(k), ALLOCATABLE,INTENT(OUT) :: yarray(:)

   outfile = 'initialized.sac'

   delta=rundef; depmin=rundef; depmax=rundef; scale=rundef; odelta=rundef;
   b=rundef; e=rundef; o=rundef; a=rundef; internal1=rundef; t0=rundef;
   t1=rundef; t2=rundef; t3=rundef; t4=rundef; t5=rundef; t6=rundef;
   t7=rundef; t8=rundef; t9=rundef; f=rundef; resp0=rundef; resp1=rundef;
   resp2=rundef; resp3=rundef; resp4=rundef; resp5=rundef; resp6=rundef;
   resp7=rundef; resp8=rundef; resp9=rundef; stla=rundef; stlo=rundef;
   stel=rundef; stdp=rundef; evla=rundef; evlo=rundef; evel=rundef;
   evdp=rundef; mag=rundef; user0=rundef; user1=rundef; user2=rundef;
   user3=rundef; user4=rundef; user5=rundef; user6=rundef; user7=rundef;
   user8=rundef; user9=rundef; dist=rundef; az=rundef; baz=rundef;
   gcarc=rundef; internal2=rundef; internal3=rundef; depmen=rundef;
   cmpaz=rundef; cmpinc=rundef; xminimum=rundef; xmaximum=rundef;
   yminimum=rundef; ymaximum=rundef; unused1=rundef; unused2=rundef;
   unused3=rundef; unused4=rundef; unused5=rundef; unused6=rundef;
   unused7=rundef;

   nzyear=iundef; nzjday=iundef; nzhour=iundef; nzmin=iundef; 
   nzsec=iundef; nzmsec = iundef;
   nvhdr=6; norid=iundef; nevid=iundef; npts=iundef; internal4=iundef;
   nwfid=iundef; nxsize=iundef; nysize=iundef; unused8=iundef
   iftype=1; idep=iundef; iztype=iundef; unused9=iundef; iinst=iundef;
   istreg=iundef; ievreg=iundef; ievtyp=iundef; iqual=iundef; isynth=iundef;
   imagtyp=iundef; imagsrc=iundef; unused10=iundef; unused11=iundef;
   unused12=iundef; unused13=iundef; unused14=iundef; unused15=iundef;
   unused16=iundef; unused17=iundef; 
   
   leven=1; lpspol=0; lovrok=1; lcalda=1; unused18=0;

   kevnm=kundef; kstnm=kundef; khole=kundef; ko=kundef;
   ka=kundef; kt0=kundef; kt1=kundef; kt2=kundef; kt3=kundef;
   kt4=kundef; kt5=kundef; kt6=kundef; kt7=kundef; kt8=kundef
   kt9=kundef; kf=kundef; kuser0=kundef; kuser1=kundef; kuser2=kundef
   kcmpnm=kundef; knetwk=kundef; kdatrd=kundef; kinst=kundef

   ALLOCATE(yarray(1))
   yarray = 0.0;


  END SUBROUTINE initsac

END MODULE sac_i_o
