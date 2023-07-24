PWD_fold=`pwd`
DATA_fold=${PWD_fold}/data/
StaCode=$1
FIG_fold=${PWD_fold}/figs/
TMP_fold=${PWD_fold}/tmp/

cd $TMP_fold
psfile=$FIG_fold/$StaCode.Fig07.rmse.ps
xyz_b=$TMP_fold/grd.b.xyz
xyz_a=$TMP_fold/grd.a.xyz
grd_b=$TMP_fold/b.grd
grd_a=$TMP_fold/a.grd
rmse_b=$TMP_fold/rmse.b.xy
rmse_a=$TMP_fold/rmse.a.xy
reft=$TMP_fold/reft.xy
cat $DATA_fold$StaCode/${StaCode}SRF_before.xyz |awk '{printf("%.1f %d %f\n",$1,$2*1,$3)}'>$xyz_b
cat $DATA_fold$StaCode/${StaCode}SRF_after.xyz |awk '{printf("%.1f %d %f\n",$1,$2*1,$3)}'>$xyz_a
cat $DATA_fold/$StaCode/${StaCode}rmse_before.xy |awk '{printf("%d %f\n",$1,$2)}'>$rmse_b
cat $DATA_fold/$StaCode/${StaCode}rmse_after.xy |awk '{printf("%d %f\n",$1,$2)}'>$rmse_a
cat $DATA_fold/$StaCode/${StaCode}.ref_sum.xy >$reft
RangeB=`cat $xyz_b|gmt gmtinfo -C|awk '{print $1"/"$2-45"/"$3"/"$4}'`
RangeA=`cat $xyz_a|gmt gmtinfo -C|awk '{print $1"/"$2-45"/"$3"/"$4}'`
RangeNum=`cat $rmse_a|gmt gmtinfo -C|awk '{print $1"/"$2}'`
NumAll=`cat $xyz_a|gmt gmtinfo -C|awk '{print $4}'`
gmt xyz2grd $xyz_b -R$RangeB -G$grd_b -I0.1/1
gmt xyz2grd $xyz_a -R$RangeA -G$grd_a -I0.1/1

gmt gmtset PS_MEDIA A2
gmt psbasemap -R$RangeB -JX3i/6.5i -Bx10f10+l"Time (s)" -By10f5 -BWnS -P -K>$psfile
gmt makecpt -Cpolar -T-0.4/0.4 -D> tmp.cpt
gmt grdimage  $grd_b  -J -R -B -E1000 -Ctmp.cpt -O -K>> $psfile
gmt pstext -R -J -F+f18p,1,black+jRT -Ggray -B -N -K -O>>$psfile<< EOF
`echo $NumAll |awk '{printf("48 %f RAW",$1-5)}'`
EOF
gmt psbasemap -R-0.2/1/$RangeNum -JX1i/6.5i -Bx0.4f0.2+l"RMSE" -By10f5 -BewnS -X3i -O -K>>$psfile

cat $rmse_b|awk '{print $2,$1}'|gmt psxy -J -R -W1p,black -O -K >> $psfile
cat $rmse_a|awk '{print $2,$1}'|gmt psxy -J -R -W1p,red -O -K>> $psfile
gmt pstext -R -J -F+f12p,1,red+jRT -B -N -K -O>>$psfile<< EOF
`echo $NumAll|awk '{printf("0.9 %f Resorted",$1-$1*0.05)}'`
EOF
gmt pstext -R -J -F+f12p,1,black+jRT -B -N -K -O>>$psfile<< EOF
`echo $NumAll|awk '{printf("0.9 %f RAW",$1*0.10)}'`
EOF

gmt psbasemap -R$RangeA -JX3i/6.5i -Bx10f10+l"Time (s)" -By10f5 -BneS -X1i -P -K -O >>$psfile
gmt grdimage  $grd_a  -J -R -B -E1000 -Ctmp.cpt -O -K>> $psfile
gmt pstext -R -J -F+f18p,1,red+jRT -Ggray -B -N -K -O>>$psfile<< EOF
`echo $NumAll |awk '{printf("48 %f Resorted",$1-5)}'`
EOF
gmt psbasemap -R0/9/-0.05/0.12 -JX2i/2i -Bx1f1+l"Time (s)" -By0.02f0.01 -BeWNs -X3.8i -Y4.5i -O -K>>$psfile
cat $reft|awk '{print $1,$2}'|gmt psxy -J -R -W0.5p,gray -O -K>> $psfile
cat $reft|awk '{print $1,$3}'|gmt psxy -J -R -W0.5p,black,- -O -K>> $psfile
cat $reft|awk '{print $1,$4}'|gmt psxy -J -R -W0.5p,black,..- -O -K>> $psfile
cat $reft|awk '{print $1,$5}'|gmt psxy -J -R -W0.5p,blue,- -O -K>> $psfile
cat $reft|awk '{print $1,$6}'|gmt psxy -J -R -W0.5p,brown,-  -O -K>> $psfile
cat $reft|awk '{print $1,$7}'|gmt psxy -J -R -W0.5p,blue,..- -O -K>> $psfile
cat $reft|awk '{print $1,$8}'|gmt psxy -J -R -W1p,red -O -K>> $psfile
cat $reft|awk '{print $1,$9}'|gmt psxy -J -R -W0.5p,brown,..- -O -K>> $psfile
cat $reft|awk '{print $1,$10}'|gmt psxy -J -R -W0.5p,darkgreen,- -O -K>> $psfile
cat $reft|awk '{print $1,$11}'|gmt psxy -J -R -W0.5p,darkgreen,..- -O -K>> $psfile
cat $reft|awk '{print $1,$12}'|gmt psxy -J -R -W0.5p,cyan,- -O -K>> $psfile
cat $reft|awk '{print $1,$13}'|gmt psxy -J -R -W0.5p,cyan,..- -O -K>> $psfile
cat $reft|awk '{print $1,$14}'|gmt psxy -J -R -W0.5p,pink,- -O -K>> $psfile
gmt pstext -R -J -F+f+a+j -Ggray -B -N -K -O>>$psfile<< EOF
0.5 0.11 18 0 0 LT Moho
EOF

gmt psbasemap -R10/25/-0.12/0.05 -JX2i/2i -Bx2f1+l"Time (s)" -By0.02f0.01 -BeWNs -X2.8i -O -K>>$psfile
cat $reft|awk '{print $1,$2}'|gmt psxy -J -R -W0.5p,gray -O -K>> $psfile
cat $reft|awk '{print $1,$3}'|gmt psxy -J -R -W0.5p,black,- -O -K>> $psfile
cat $reft|awk '{print $1,$4}'|gmt psxy -J -R -W0.5p,black,..- -O -K>> $psfile
cat $reft|awk '{print $1,$5}'|gmt psxy -J -R -W0.5p,blue,- -O -K>> $psfile
cat $reft|awk '{print $1,$6}'|gmt psxy -J -R -W0.5p,brown,-  -O -K>> $psfile
cat $reft|awk '{print $1,$7}'|gmt psxy -J -R -W0.5p,blue,..- -O -K>> $psfile
cat $reft|awk '{print $1,$8}'|gmt psxy -J -R -W1p,red -O -K>> $psfile
cat $reft|awk '{print $1,$9}'|gmt psxy -J -R -W0.5p,brown,..- -O -K>> $psfile
cat $reft|awk '{print $1,$10}'|gmt psxy -J -R -W0.5p,darkgreen,- -O -K>> $psfile
cat $reft|awk '{print $1,$11}'|gmt psxy -J -R -W0.5p,darkgreen,..- -O -K>> $psfile
cat $reft|awk '{print $1,$12}'|gmt psxy -J -R -W0.5p,cyan,- -O -K>> $psfile
cat $reft|awk '{print $1,$13}'|gmt psxy -J -R -W0.5p,cyan,..- -O -K>> $psfile
cat $reft|awk '{print $1,$14}'|gmt psxy -J -R -W0.5p,pink,- -O -K>> $psfile
gmt pstext -R -J -F+f+a+j -Ggray -B -N -K -O>>$psfile<< EOF
10.5 0.04 18 0 0 LT LAB
EOF
gmt pslegend -R -J -F+p0p -D+w2i/2i+l0.8i -X2.5i -O -K>> $psfile << EOF
S 0.8c - 0.4i - 0.5p,gray 1.5c 100%(rmse ref)
S 0.8c - 0.4i - 0.5p,black,- 1.5c 5% 
S 0.8c - 0.4i - 0.5p,black,..- 1.5c 10%
S 0.8c - 0.4i - 0.5p,blue,- 1.5c 15%
S 0.8c - 0.4i - 0.5p,brown,- 1.5c 20%
S 0.8c - 0.4i - 0.5p,blue,..- 1.5c 25%
S 0.8c - 0.4i - 1p,red 1.5c 30% 
S 0.8c - 0.4i - 0.5p,brown,..- 1.5c 40%
S 0.8c - 0.4i - 0.5p,darkgreen,- 1.5c 50%
S 0.8c - 0.4i - 0.5p,darkgreen,..- 1.5c 60%
S 0.8c - 0.4i - 0.5p,cyan,- 1.5c 70%
S 0.8c - 0.4i - 0.5p,cyan,..- 1.5c 80%
S 0.8c - 0.4i - 0.5p,pink,- 1.5c 90%
EOF

rm -rf sac_tmp
mkdir sac_tmp
pltlst=pssac.plt.lst
cp $DATA_fold$StaCode/*rmse30 ./sac_tmp/
saclst user1 f ./sac_tmp/*rmse30|sort -k2|awk '{print $1,-5,NR}'> $pltlst
NumSRF=`cat $pltlst|wc|awk '{print $1+2}'`
Num30=`cat $pltlst|wc|awk '{print $1}'`
gmt psbasemap -R-5/50/-2/$NumSRF -JX4.8i/3.8i -Bx10f10+l"Time (s)" -By10f5 -BeWnS -X-5.3i -Y-4.5i -O -K>>$psfile
#saclst user1 f ./sac_tmp/*rmse30|sort -k2|awk '{print $2,NR}'> pssac.plt.lst
cat pssac.plt.lst|gmt pssac -R -J -M0.8 -W0.5p,black -Gn+gdeepskyblue1 -Gp+gfirebrick1 -B -O -K>> $psfile

BIN=$PWD_fold/../Section001/bin
${BIN}/sacStack -N -O$StaCode.stack2021.sac<<eof
`saclst user1 f ./sac_tmp/*rmse30|sort -k2|awk '{print $1,0,1}' `
eof
gmt psbasemap -R-5/50/-1/1 -JX4.8i/0.5i -Bx20f10g10 -By10f5 -Bnsew -O -K -Y3.8i >>$psfile
echo $StaCode.stack2021.sac 0|gmt pssac -R -J -M0.8 -W0.5p,black -Gn+gdeepskyblue1 -Gp+gfirebrick1 -B -O -K>> $psfile

gmt pstext -R -J -F+f+a+j -B -N -K -O>>$psfile<< EOF
47.5 0.9 12 0 0 RT Sum rmse 30% 
55 1 16 0 0 LT Station: $StaCode
55 0 16 0 0 LT RMSE 30%: ($Num30/$NumAll)
EOF


gmt psbasemap -R0/120/-2/$NumSRF -JX0.8i/3.8i -Bx60f30g30+l"GCARC" -By10f5 -BwS -O -K -X5.3i -Y-3.8i >>$psfile
saclst user1 gcarc baz f ./sac_tmp/*rmse30| sort -k2|awk '{print $3,NR}'|gmt psxy -R -J -W0.25p,black -Gred -B -St5p -O -K>> $psfile
gmt psbasemap -R0/360/-2/$NumSRF -JX0.8i/3.8i -Bx180f90g90+l"BAZ" -By10f5 -BSw -O -K -X1.2i>>$psfile
saclst user1 gcarc baz f ./sac_tmp/*rmse30| sort -k2|awk '{print $4,NR}'|gmt psxy -R -J -W0.25p,black -Gblue -B -St5p -O >> $psfile

# echo|gmt psxy -J -R -O >> $psfile
cd $FIG_fold
cp $psfile $psfile.bk.ps
gmt psconvert -A $psfile