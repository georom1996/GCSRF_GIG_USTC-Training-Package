export SAC_DISPLAY_COPYRIGHT=0
pwd_path=$(pwd)
wk_fd=${pwd_path}/../example01_kumar
BIN=${pwd_path}/../bin
fig_fd=${pwd_path}/./figs/
sacdatapath=$wk_fd/data
stainfofile=$wk_fd/stationinfo.lst
cd $sacdatapath
# XR.ST09.01.BHZ.M.1998.197.121411.SAC
zcom_filename=XR.ST09.01.BHZ.M.1998.197.121411.SAC
ncom_filename=$(echo $zcom_filename | awk 'split($0,aa,".BHZ."){print aa[1]".BHN."aa[2]}')
ecom_filename=$(echo $zcom_filename | awk 'split($0,aa,".BHZ."){print aa[1]".BHE."aa[2]}')
staname=$(echo $zcom_filename | awk 'split($0,aa,"."){print aa[1]"."aa[2]}')
stanm=$(echo $zcom_filename | awk 'split($0,aa,"."){print aa[2]}')
evttime=$(echo $zcom_filename | awk 'split($0,aa,"."){print aa[6]"."aa[7]"."aa[8]}')
evtday=$(echo $zcom_filename | awk 'split($0,aa,"."){print aa[7]}')
setname=$staname.$evttime
if [ ! -d $wk_fd/calculate_test ];then
    mkdir $wk_fd/calculate_test
else
    echo $wk_fd/calculate_test exist
fi    
cd $wk_fd/calculate_test
cp $sacdatapath/$zcom_filename ./
cp $sacdatapath/$ncom_filename ./
cp $sacdatapath/$ecom_filename ./
stlainfo=$(cat $stainfofile | grep $stanm | awk 'split($0,aa,","){print aa[2]}')
stloinfo=$(cat $stainfofile | grep $stanm | awk 'split($0,aa,","){print aa[3]}')
stelinfo=$(cat $stainfofile | grep $stanm | awk 'split($0,aa,","){print aa[4]}')
evlainfo=$(echo $evtday | awk '{if($1=="028"){print 52.886}else if($1=="197"){print -11.04}}')
evloinfo=$(echo $evtday | awk '{if($1=="028"){print -169.123}else if($1=="197"){print 166.16}}')
evdpinfo=$(echo $evtday | awk '{if($1=="028"){print 67.2}else if($1=="197"){print 110.2}}')
gmtinfo=$(echo $evtday | awk '{if($1=="028"){printf("%4d %03d %02d %02d %02d %03d\n",1999,28,8,10,5,420)}else if($1=="197"){printf("%4d %03d %02d %02d %02d %03d\n",1998,197,11,56,36,420)}}')
origininfo=$(echo $evtday | awk '{if($1=="028"){printf("%s\n","JAN-28-(28)-1999_08:10:05.420")}else if($1=="197"){printf("%s\n","JUL-16-(197)-1998_11:56:36.420")}}')
maginfo=$(echo $evtday | awk '{if($1=="028"){printf("%s\n","6.6")}else if($1=="197"){printf("%s\n","7.0")}}')
sac <<eof
r $zcom_filename $ncom_filename $ecom_filename
ch stla $stlainfo
wh
ch stlo $stloinfo
wh
ch stel $stelinfo
wh
ch evla $evlainfo
wh
ch evlo $evloinfo
wh
ch evdp $evdpinfo
wh
w over
q
eof
gcarcinfo=$(saclst gcarc f $zcom_filename | awk '{print $2}')
bazinfo=$(saclst baz f $zcom_filename | awk '{print $2}')
sac <<eof
r $zcom_filename
ch o gmt $gmtinfo
ch allt (0-&1,o&) iztype IO
wh
w over
r $ncom_filename
ch o gmt $gmtinfo
ch allt (0-&1,o&) iztype IO
wh
w over
r $ecom_filename
ch o gmt $gmtinfo
ch allt (0-&1,o&) iztype IO
wh
w over
q
eof
cp $wk_fd/iasp91.dat ./
taup_setsac -mod iasp91 -ph S-9 -evdpkm $zcom_filename
taup_setsac -mod iasp91 -ph S-9 -evdpkm $ncom_filename
taup_setsac -mod iasp91 -ph S-9 -evdpkm $ecom_filename
sac <<eof
cut t9 -200 200
r $zcom_filename $ncom_filename $ecom_filename
w $setname.z.sac $setname.n.sac $setname.e.sac
cut off
r $setname.z.sac $setname.n.sac $setname.e.sac
rtrend
rmean
lp c 0.25
interpolate delta 0.1
w over
q
eof
for inci in 00 04 08 12 16 20 24 28 32 36 40 44 48 52 56 60; do
    $BIN/SX_Deconv -Z $setname.z.sac -N $setname.n.sac -E $setname.e.sac -t1 190 -t2 280 -reg 1 -inc $inci -D
    mv fil.sac fil.$setname.$inci.sac
    mv l-com.sac l-com.$setname.$inci.sac
    mv q-com.sac q-com.$setname.$inci.sac
    mv t-com.sac t-com.$setname.$inci.sac
    mv rf.nsl.sac srf.l.$setname.$inci.sac
    mv rf.nsq.sac srf.q.$setname.$inci.sac
    mv rf.nst.sac srf.t.$setname.$inci.sac
    sac <<eof
cut -100 100
r srf.l.$setname.$inci.sac
w srf.l.$setname.$inci.sac.100.cut
cut off
r srf.l.$setname.$inci.sac.100.cut
reverse
w over
q
eof
    $BIN/spmove.joe <<eof
srf.l.$setname.$inci.sac.100.cut
srf.l.$setname.$inci.sac.100.cut.rm
1
eof
    sac <<eof
r srf.l.$setname.$inci.sac.100.cut.rm
reverse
w over
q
eof
done
sac <<eof
cut -30 30
r srf.l.$setname.*.sac.100.cut.rm
w append .cut
cut off
q
eof
ls srf.l.$setname.*.sac.100.cut.rm.cut | awk '{print $1,(NR-1)*4,'-30',"1p,black"}' >lsrf.$setname.lst.new

gmt pssac lsrf.$setname.lst.new -JX6i/2.5i -R-3.9/63.9/-30/0 -Bxa4f4+l"corrected - $staname - Incidence Angle (@~\260@~)" -Bya10f5g5+l"Time (s)" -BWSne -M1 -Fr -C-30/0 -Gn+g0/0/0+t-30/0+z0 -Gp+ggray+t-30/0 -Q -P -K >$fig_fd/corrected-$setname-SRF.l.ps
gmt pssac -JX4i/2i -R-300/200/0.5/3.5 -Bxa100f50+l"Time (s)" -BWS -M1 -Fr -T+t9 -W1p,black -C-300/200 -P -Y3.2i -K -O >>$fig_fd/corrected-$setname-SRF.l.ps <<eof
$zcom_filename -300 1
$ncom_filename -300 2
$ecom_filename -300 3
eof
gmt psxy -J -R -P -O -K -W1p,red,- >>$fig_fd/corrected-$setname-SRF.l.ps <<eof
0 0.5
0 3.5
eof
gmt psxy -J -R -P -O -K -W1p,gray >>$fig_fd/corrected-$setname-SRF.l.ps <<eof
-10 0.5
-10 3.5
eof
gmt psxy -J -R -P -O -K -W1p,gray >>$fig_fd/corrected-$setname-SRF.l.ps <<eof
80 0.5
80 3.5
eof
gmt pstext -J -R -P -O -K >>$fig_fd/corrected-$setname-SRF.l.ps <<eof
-50 3.3 S-onset
eof

gmt pstext -JX0.5i/2i -R-1/1/0.5/3.5 -P -O -K -X-0.4i >>$fig_fd/corrected-$setname-SRF.l.ps <<eof
0 1 Z
0 2 N
0 3 E
eof
gmt pstext -JX3i/2i -R0/10/0.5/3.5 -P -O -X4.45i -Y-0.25i -F+f10p,1,black+a0+jLT >>$fig_fd/corrected-$setname-SRF.l.ps <<eof
0 3.4 Station: $staname
0 3.0 Origin: 
0 2.6 $origininfo
0 2.2 Magnitude: $maginfo
0 1.8 Distance: $gcarcinfo
0 1.4 BakAzimuth: $bazinfo
eof
gmt psconvert -A $fig_fd/corrected-$setname-SRF.l.ps
