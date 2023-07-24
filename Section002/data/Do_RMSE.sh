export SAC_DISPLAY_COPYRIGHT=0
workpath=$(pwd)
cd $workpath
for stnm in $(cat fdd.lst); do
	cd $stnm
	sac <<eof
r *rm
reverse
w append .rv
cut -5 95
r *.rm.rv
w append .cut
cut off
q
eof
# 	for sacname in $(ls *rm.rv.cut); do
# 		deg=$(saclst gcarc f $sacname | awk '{print $2}')
# 		evdp=$(saclst evdp f $sacname | awk '{print $2}')
# 		rayp=$(taup_time -mod iasp91 -h $evdp -deg $deg -ph S | awk 'NR==6{print $5/111}')
# 		sac <<eof
# r $sacname
# ch user0 $rayp
# wh
# q
# eof
# 	done
	cd $workpath
done

cd $workpath
for stnm in $(cat fdd.lst); do
	cd $stnm
	for sacname in $(ls *rm.rv.cut); do
		newname=$(echo $sacname | awk '{split($1,aa,".sac.");print aa[1]".srf"}')
		cp $sacname $newname
	done
	cd $workpath
done
echo "octave start >>>"
octave --no-line-editing --silent --no-gui Silent_FiltRMSE.m
echo "octave done <<<"
cd $workpath
for stnm in $(cat fdd.lst); do
	cd $stnm
	bash ./$stnm.rmse10.lst
	bash ./$stnm.rmse20.lst
	bash ./$stnm.rmse30.lst
	cd $workpath
	echo $stnm RMSE Done~
done
