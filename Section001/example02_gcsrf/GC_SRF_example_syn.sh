export SAC_DISPLAY_COPYRIGHT=0
pwd_path=$(pwd)
wk_fd=${pwd_path}/../example02_gcsrf
syn_data_fd=${pwd_path}/../example02_gcsrf/data/syn_SM_d75/
BIN=${pwd_path}/../bin

setname=2021.000.00.00.00.0000.RF.d75.00
nfilename=$setname.n.sac
efilename=$setname.e.sac
zfilename=$setname.z.sac

mkdir -p $wk_fd/calculate_test/syn_SM_d75
cd $wk_fd/calculate_test/syn_SM_d75
cp $syn_data_fd/$nfilename ./
cp $syn_data_fd/$efilename ./
cp $syn_data_fd/$zfilename ./
cp $wk_fd/iasp91.dat ./

S_on=$(saclst t9 f $zfilename | awk '{print $2}')
ts_b=$(echo $S_on | awk '{print $1-5}')
ts_e=$(echo $S_on | awk '{print $1+30}')
tn_b=$(echo $S_on | awk '{print $1-60}')
tn_e=$(echo $S_on | awk '{print $1-30}')
snr_n=$($BIN/sacsnr $nfilename -s $ts_b $ts_e -n $tn_b $tn_e | awk 'END{print $5}')
snr_e=$($BIN/sacsnr $efilename -s $ts_b $ts_e -n $tn_b $tn_e | awk 'END{print $5}')
snr_z=$($BIN/sacsnr $zfilename -s $ts_b $ts_e -n $tn_b $tn_e | awk 'END{print $5}')
SNR=$(echo $snr_z $snr_n $snr_e | awk '{print ($2+$3)/2}')
SNR_flag=$(echo $snr_z $snr_n $snr_e | awk '{if(($2+$3)/2>=5){print 1}else{print 0}}')
echo "Signal Noise Ratio: z n e avr_ne"
echo $snr_z $snr_n $snr_e $SNR
echo "Signal Noise Ratio: z n e avr_ne" >SNR_REC.DAT
echo $snr_z $snr_n $snr_e $SNR >>SNR_REC.DAT
echo $SNR_flag >>SNR_REC.DAT

if [ $SNR_flag -eq 1 ]; then
	{

		timebegin=190
		$BIN/GC_SRF_par -Z $setname.z.sac -N $setname.n.sac -E $setname.e.sac -t1 $timebegin -reg 3

		sac <<eof
cut -100 100
r rf*nsl.sac
w append .100.cut
cut off
r *.100.cut
reverse
w append rv
q
eof

		for srfnm in $(ls *100.cutrv); do
			$BIN/spmove.joe <<eof
$srfnm
$srfnm.rm
1
eof
		done

		sac <<eof
r *100.cutrv.rm
reverse
w append .rv
cut -30 30
r *100.cutrv.rm.rv
w append cut
cut off
q
eof
		#Post CCCCCCCCCCCCCCCCCCCCC
		$BIN/GC_SRF_post
		#Post CCCCCCCCCCCCCCCCCCCCC
		gcarcinfo=$(saclst gcarc f $zfilename | awk '{print $2*1}')
		evdpinfo=$(saclst evdp f $zfilename | awk '{print $2*1}')
		#taup_time -mod iasp91 -ph S -deg $gcarcinfo -h $evdpinfo|awk 'NR==6{print $7*1}'
		inci_syn=$(taup_time -mod iasp91 -ph S -deg $gcarcinfo -h $evdpinfo | awk 'NR==6{print $7*1}')
		echo "See: Gcarc Evdp Inci_syn"
		echo $gcarcinfo $evdpinfo $inci_syn
		echo $gcarcinfo $evdpinfo $inci_syn >Syn.dat
	}
fi
