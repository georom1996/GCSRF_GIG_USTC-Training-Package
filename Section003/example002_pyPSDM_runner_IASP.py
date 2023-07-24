# Step1 : data preparation ::
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
import obspy
import shutil
import pandas as pd
import numpy as np
import subprocess
from subprocess import PIPE, Popen
from obspy.taup import TauPyModel
import os
from scipy import interpolate
# from obspy.taup.velocity_model import VelocityModel
# from obspy.taup.taup_create import build_taup_model
from obspy import read
import itertools
import os
import datetime
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")
cwdfd = os.getcwd()
Ori_optpath = cwdfd+"/data_gcsrf/Opt_Papermods/"
stalst_file = cwdfd+"/psdm/PSDM_qseis2022_IASP/SP.lat.lon.lst"
Mod_path = cwdfd+"/mod/"
PSDM_DataPath = cwdfd+"/psdm/PSDM_qseis2022_IASP/data/"
PSDM_pro_path = cwdfd+"/psdm/PSDM_qseis2022_IASP/"
subname = "QseisMig"

# srfname='CCPUSE.srf.l.2016.022.16.44.49.4000.SC.SMK.00.36-60.sac.rv.cut'


def calc_tsp_mout(modelname, H_depth, rayp_ref):
    datamod = pd.read_csv(modelname, sep='\s+', header=None,
                          names=['depth', 'vp', 'vs', 'rho'], skiprows=lambda x: x in [0, 1])
    Depths = datamod.depth
    Vs = datamod.vs
    Vp = datamod.vp
    Rho = datamod.rho
    dz = np.append(0, np.diff(Depths))
    RR = 6371 - Depths
    idx_top = []
    for ii in range(len(Depths)):
        if Depths[ii] < H_depth:
            idx_top.append([ii, Depths[ii], Vp[ii], Vs[ii]])
    idx_top.append([len(idx_top), H_depth, idx_top[-1][2], idx_top[-1][3]])
    TS_Sp_ref = 0
    for ll in range(len(idx_top)-1):
        TS_Sp_ref += (idx_top[ll+1][1]-idx_top[ll][1])*(np.sqrt(np.abs(1/(idx_top[ll+1][3])
                                                                       ** 2-rayp_ref**2))-np.sqrt(np.abs(1/(idx_top[ll+1][2])**2-rayp_ref**2)))
    return TS_Sp_ref


def moveout_joe(srfname):
    outname = srfname+'.mout'
    rayp_ref = 0.06
    depth_begin = 1
    depth_end = 1200
    velocity_model = Mod_path+"iasp91.tvel"
    st = read(srfname)
    # st[0].plot()
    tr = st[0]
    trout = tr.copy()
    timeseries = np.linspace(tr.stats.sac.b, tr.stats.sac.e, tr.stats.sac.npts)
    ampseries = tr.data
    # filename = "./iasp91.tvel"
    # rayp = model.get_ray_paths(source_depth_in_km=0,distance_in_degree=tr.stats.sac.gcarc, phase_list=['S'])
    # S_rayp=rayp[0].ray_param*np.pi/180/111.591
    S_rayp = tr.stats.sac.user0
    ts_old = []
    ts_new = []
    for dp in np.linspace(depth_begin, depth_end, depth_end):
        ts_old.append(calc_tsp_mout(velocity_model, dp, S_rayp))
        ts_new.append(calc_tsp_mout(velocity_model, dp, rayp_ref))
    # print(np.array(ts_old)-np.array(ts_new))
    timeseries_new = interpolate.interp1d(
        ts_old, ts_new, bounds_error=False, fill_value='extrapolate')(timeseries)
    # ampseries_new=interpolate.interp1d(timeseries_new, ampseries,bounds_error=False, fill_value='extrapolate')(timeseries_new)
    timeseries_new_int = np.linspace(
        tr.stats.sac.b, tr.stats.sac.e, tr.stats.sac.npts)
    ampseries_new_int = interpolate.interp1d(
        timeseries_new, ampseries, bounds_error=False, fill_value='extrapolate')(timeseries_new_int)
    trout.stats.channel = 'mout.joe'
    trout.stats.sampling_rate = tr.stats.sampling_rate
    trout.stats.delta = tr.stats.delta
    trout.stats.sac.delta = tr.stats.sac.delta
    trout.npts = len(ampseries_new_int)
    trout.stats.sac.b = tr.stats.sac.b
    trout.stats.sac.e = tr.stats.sac.e
    trout.data = ampseries_new_int
    # trout.plot()
    trout.write(outname, format='SAC')
    # print(trout.stats)


# cmd
first_run = False
stalstdata = pd.read_csv(stalst_file, delim_whitespace=True, header=None, names=[
                         'stnm', 'stlo', 'stla', 'stel', 'moho', 'lab'])
os.chdir(PSDM_pro_path)
os.getcwd()
sub_name = "QseisMig"
model = TauPyModel(model="iasp91")
if (first_run):
    if not os.path.isdir(PSDM_DataPath+subname):
        os.mkdir(PSDM_DataPath+subname)

    Sta_namelst = []
    count_rm = 0
    for names in sorted(os.listdir(Ori_optpath)):
        if (names.endswith("nsl.sac") and (names.split("_")[1] in list(stalstdata['stnm']))):
            Sta_namelst.append(names.split("_")[1])
            if (names.split(".")[11] == "00" or names.split(".")[12] == "000" or names.split(".")[11] == "60" or float(names.split(".")[12]) >= 100):
                # shutil.move(Ori_optpath+names,Ori_optpath+names+".removed")
                count_rm += 1
    Uniq_lst = sorted(list(set(Sta_namelst)))
    print(Uniq_lst)
    countmax = 0
    for staname in Uniq_lst:

        if not os.path.isdir(PSDM_DataPath+subname+"/"+staname):
            os.mkdir(PSDM_DataPath+subname+"/"+staname)
        for sac_names in sorted(os.listdir(Ori_optpath)):
            if (sac_names.endswith("nsl.sac") and (staname in sac_names)):
                shutil.copy(Ori_optpath+"/"+sac_names,
                            PSDM_DataPath+subname+"/"+staname+"/"+sac_names)
        os.chdir(PSDM_DataPath+subname+"/"+staname)
        os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
        sac_cmd = "cut -100 100\n"
        sac_cmd += "r *nsl.sac\n"
        sac_cmd += "reverse\n"
        sac_cmd += "w append .cut100rv\n"
        sac_cmd += "cut off\n"
        sac_cmd += "cut -5 95\n"
        sac_cmd += "r *nsl.sac.cut100rv\n"
        sac_cmd += "w over\n"
        sac_cmd += "cut off\n"
        sac_cmd += "q\n"
        subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(
            sac_cmd.encode())
        for names in sorted(os.listdir(PSDM_DataPath+subname+"/"+staname)):
            if (names.endswith("cut100rv")):
                st0 = obspy.read(PSDM_DataPath+subname+"/" +
                                 staname+"/"+names, headonly=True)
                if (st0[0].stats.sac.depmax > 1):
                    shutil.move(PSDM_DataPath+subname+"/"+staname+"/"+names,
                                PSDM_DataPath+subname+"/"+staname+"/"+names+".removed")
                else:
                    gcarc_0 = st0[0].stats.sac.gcarc
                    evdp_0 = st0[0].stats.sac.evdp
                    S_arrivals = model.get_travel_times(
                        source_depth_in_km=evdp_0, distance_in_degree=gcarc_0, phase_list=["S"])
                    ray_names = S_arrivals[0].ray_param*2*np.pi/360/111.195
                    os.chdir(PSDM_DataPath+subname+"/"+staname+"/")
                    # sac_rayp_cmd = "cd "+PSDM_DataPath+subname+"/"+staname+"/\n"
                    sac_rayp_cmd = "r "+names+"\n"
                    sac_rayp_cmd += "ch user0 "+str(ray_names)+"\n"
                    sac_rayp_cmd += "wh\n"
                    sac_rayp_cmd += "q\n"
                    # print(sac_rayp_cmd)
                    subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(
                        sac_rayp_cmd.encode())
    for staname in Uniq_lst:
        os.chdir(PSDM_DataPath+subname+"/"+staname)
        # print(PSDM_DataPath+subname+"/"+staname)
        os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
        sac_cmd = "r *.cut100rv\n"
        sac_cmd += "rmean\n"
        sac_cmd += "rtr\n"
        sac_cmd += "w over\n"
        sac_cmd += "q\n"
        subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(
            sac_cmd.encode())

    for staname in Uniq_lst:
        os.chdir(PSDM_DataPath+subname+"/"+staname)
        srfname_lst = []
        for names in sorted(os.listdir(PSDM_DataPath+subname+"/"+staname)):
            if (names.endswith("sac.cut100rv")):
                srfname_lst.append(names)
        start_t1 = datetime.datetime.now()
        # print(srfname_lst)
        num_cores = int(mp.cpu_count())
        pool = mp.Pool(num_cores)
        final = pool.map(moveout_joe, srfname_lst)
        end_t1 = datetime.datetime.now()
        elapsed_sec1 = (end_t1 - start_t1).total_seconds()
        print("Moveout Corrected at %s Cost: %.4f s" % (staname, elapsed_sec1))


os.chdir(PSDM_DataPath+subname)
proc_GetLst = Popen(
    'bash '+PSDM_DataPath+"get_datalist.qseis.mout.sh",
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
outinfo_GetLst, errinfo_GetLst = proc_GetLst.communicate()

######## Parameters ALL Here ########

##################################################################################################
# [m660q Parameters Part Begin]
ref_model = 'iasp.joe'
Ps_ray = 'mray_joe.Ps.dat'
Sp_ray = 'mray_joe.Sp.dat'
m6690q_Ps_out = 'm660q_Pcs.out'
m6690q_Sp_out = 'm660q_Scp.out'
iflat = 0
itype_Ps = 1
itype_Sp = -1
# Parameters I concerned about
# build mray_joe.Ps.dat and mray_joe.Sp.dat
LayerCount = 51
# [m660q Parameters Part End]

##################################################################################################
# [Pierc_new_n Parameters Part Begin]
ref_model = '../model/iasp.joe'
Ps_ray = 'mray_joe.Ps.dat'
Sp_ray = 'mray_joe.Sp.dat'
m660q_Ps_out = 'm660q_Pcs.out'
m660q_Sp_out = 'm660q_Scp.out'
pierc_out = 'pierc_use.dat'
center_la = 27.5
center_lo = 104
out_npts = 1251
sac_user_num_rayp = 0
event_filt_flag = 1  # 0: dist; 1: gcarc; 2: baz
event_filt_min = 70.0
event_filt_max = 80.0
rfdata_path = '../data/'
num_sub = 1
name_sub = sub_name
name_lst = 'datalist.txt'
# [Pierc_new_n Parameters Part End]

##################################################################################################
# [binr_vary_scan_n Parameters Part Begin]
Descar_la_begin = -500.0
Descar_lo_begin = 0
Descar_la_end = -500.0
Descar_lo_end = 0
Descar_step = 100
Profile_len = 1000.0
az_min = 0.0
az_max = 0.0
az_step = 15.0
bins_step = 50
trace_num_min = 2
ratio_trace = 1.0
UTM_zone = 48
out_trace_npts = 1001
out_trace_dt = 0.1
temp_folder = '../temp'
moveout_flag = 0
moveout_gcarc = 0
norm_flag = 1
stack_flag = 1
npief = 1
out_idx_flag = 0
binr_out_name = '2022_binr_out_Qseis_1s'
# [binr_vary_scan_n Parameters Part End]

##################################################################################################
# [Hdpmig.x Parameters Part Begin]
method_flag = 0
refvel_flag = 0
vscale = 1.0
freq_min = 0.03
freq_max = 0.5
freq_idx_b = 0
freq_idx_e = 40
nxmod = 21
nzmod = 1600
nx = 1024
nz = 1600
dx = 50
dz = 0.5
ntrace = 21
nt = 1001
dt = 0.1
nt0 = 2048
nt_step = 1
FD_ang = 45
smoothpt_left = 4
smoothpt_right = 4
mod_format_flag = 0
mig_model = "../model/iasp.joe"
stacked_field = "../stack/stack_"+binr_out_name+".dat"
hdpmig_out_field = sub_name+"_" + \
    str(freq_min)+"_"+str(freq_max)+"_hdpmig.joe.dat"
# [Hdpmig.x Parameters Part End]

Ps_rayfile = open((PSDM_pro_path+'m660q/'+Ps_ray), 'w+')
Sp_rayfile = open((PSDM_pro_path+'m660q/'+Sp_ray), 'w+')
# str_idx=("%2d " % LayerCount)
# str_flg='   '
# for Pcs
# build m660q_model.in for Ps
m660q_conf = open((PSDM_pro_path+'m660q/m660q_model.in'), 'w+')
print('* velocity model file', file=m660q_conf)
print("%s" % ref_model, file=m660q_conf)
print('* ray file', file=m660q_conf)
print("%s" % Ps_ray, file=m660q_conf)
print('* output file', file=m660q_conf)
print("%s" % m6690q_Ps_out, file=m660q_conf)
print('* iflat, itype (= 0: free-surface refl.; else: conversion (>0: Ps; <0: Sp))', file=m660q_conf)
print("%d     %d" % (iflat, itype_Ps), file=m660q_conf)
m660q_conf.close()
print(LayerCount, file=Ps_rayfile)
for cc_n in range(0, LayerCount):
    str_idx = ("%2d " % LayerCount)
    str_flg = '   '
    for ii in range(LayerCount, 0, -1):
        str_idx += ("%2d" % ii)
        if (ii <= cc_n):
            str_flg += ("%2d" % 3)
        else:
            str_flg += ("%2d" % 5)
    print(str_idx, file=Ps_rayfile)
    # print(str_idx)
    print(str_flg, file=Ps_rayfile)

# Run STEP 01 for Ps
os.chdir(PSDM_pro_path+'m660q')
print(">>> Step01 # Processing in folder:\n  %s\n\nM660q_Model INFO :\n" % os.getcwd())


proc010_make = Popen(
    'make clean\nmake\n',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
proc010_make.communicate()
proc010 = Popen(
    './M660q_model',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
outinfo010, errinfo010 = proc010.communicate()
# print(outinfo010.decode())
# print(errinfo010.decode())

# for Scp
# build m660q_model.in for Sp
m660q_conf = open((PSDM_pro_path+'m660q/m660q_model.in'), 'w+')
print('* velocity model file', file=m660q_conf)
print("%s" % ref_model, file=m660q_conf)
print('* ray file', file=m660q_conf)
print("%s" % Sp_ray, file=m660q_conf)
print('* output file', file=m660q_conf)
print("%s" % m6690q_Sp_out, file=m660q_conf)
print('* iflat, itype (= 0: free-surface refl.; else: conversion (>0: Ps; <0: Sp))', file=m660q_conf)
print("%d     %d" % (iflat, itype_Sp), file=m660q_conf)
m660q_conf.close()
print(LayerCount, file=Sp_rayfile)
for cc_n in range(0, LayerCount):
    str_idx = ("%2d " % LayerCount)
    str_flg = '   '
    for ii in range(LayerCount, 0, -1):

        str_idx += ("%2d" % ii)
        if (ii <= cc_n):
            str_flg += ("%2d" % 5)
        else:
            str_flg += ("%2d" % 3)
    print(str_idx, file=Sp_rayfile)
    print(str_flg, file=Sp_rayfile)

# Run STEP 01 for Sp
os.chdir(PSDM_pro_path+'m660q')
print(">>> Step01 # Processing in folder:\n  %s\n\nM660q_Model INFO :\n" % os.getcwd())
proc011_make = Popen(
    'make clean\nmake\n',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
proc011_make.communicate()
proc011 = Popen(
    './M660q_model',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
outinfo011, errinfo011 = proc011.communicate()
# print(outinfo011.decode())
# print(errinfo011.decode())


##################################################################################################
# #[Pierc_new_n Parameters Part End]

pierc_conf = open((PSDM_pro_path+'stack/pierc_new_n.in'), 'w+')
print('* output file name: iaj', file=pierc_conf)
print(("%s" % pierc_out), file=pierc_conf)
print('* the coordinate center of line: evla0,evlo0', file=pierc_conf)
print("%.2f,%.2f" % (center_la, center_lo), file=pierc_conf)
print('* output time point number: np0, irayp', file=pierc_conf)
print('%d      %d' % (out_npts, sac_user_num_rayp), file=pierc_conf)
print('* model file', file=pierc_conf)
print(("%s" % ref_model), file=pierc_conf)
print('* * ivar (0: dist; 1: gcarc; 2: baz),varmin,varmax', file=pierc_conf)
print('%d     %.1f     %.1f' %
      (event_filt_flag, event_filt_min, event_filt_max), file=pierc_conf)
print('* NW,(NWI(I),NWID(I),I=1,NW) ', file=pierc_conf)
# print('20,3 4,4 6,5 7,6 8,7 9,8 10,9 11,10 12,11 13,12 14,13 15,14 16,15 17,16 18,17 19,18 20,19 21,20 22,21 23,22 24',file=pierc_conf)

print('47,3 4,4 6,5 7,6 8,7 9,8 10,9 11,10 12,11 13,12 14,13 15,14 16,15 17,16 18,17 19,18 20,19 21,20 22,21 23,22 24,23 25,24 26,25 27,26 28,27 29,28 30,29 31,30 32,31 33,32 34,33 35,34 36,35 37,36 38,37 39,38 40,39 41,40 42,41 43,42 44,43 45,44 46,45 47,46 48,47 49,48 50,49 51', file=pierc_conf)
print('* NDW(1:5): indexs in NWI for outputting piercing points at 5 depths', file=pierc_conf)
print('6 13 16 19 36       40-, 100-, 150-, 210-, 300-km', file=pierc_conf)
print('* directory containing RFs     ', file=pierc_conf)
print("%s" % rfdata_path, file=pierc_conf)
print('* number of subdirectories    ', file=pierc_conf)
print("%d" % num_sub, file=pierc_conf)
print("%s" % name_sub, file=pierc_conf)
print("%s" % name_lst, file=pierc_conf)
pierc_conf.close()

# STEP 02
os.chdir(PSDM_pro_path+'stack')
if os.path.exists("./depth1.dat"):
    os.remove("./depth1.dat")
if os.path.exists("./depth2.dat"):
    os.remove("./depth2.dat")
if os.path.exists("./depth3.dat"):
    os.remove("./depth3.dat")
if os.path.exists("./depth4.dat"):
    os.remove("./depth4.dat")
if os.path.exists("./depth5.dat"):
    os.remove("./depth5.dat")
if os.path.exists("./sta_name.dat"):
    os.remove("./sta_name.dat")
if os.path.exists("./station.dat"):
    os.remove("./station.dat")
if os.path.exists("./events.txt"):
    os.remove("./events.txt")
if os.path.exists("./"+pierc_out):
    os.remove("./"+pierc_out)
print(">>> Step02 # Processing in folder:\n  %s\n\npierc_new_n INFO :\n" % os.getcwd())
proc02_make = Popen(
    'make clean\nmake\n',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
proc02_make.communicate()
proc02 = Popen(
    './pierc_new_n',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
outinfo02, errinfo02 = proc02.communicate()
# print(outinfo02.decode())
# print(errinfo02.decode())


##################################################################################################
# For Ps binr_vary_scan
binr_conf = open((PSDM_pro_path+'stack/binr_vary_scan_n.inp'), 'w+')
print('* begin and end coordinate of start point, point interval(km): begla0,beglo0,endla0,endlo0,dsp', file=binr_conf)
print("%.1f,%.1f,%.1f,%.1f,%.1f," % (Descar_la_begin, Descar_lo_begin,
      Descar_la_end, Descar_lo_end, Descar_step), file=binr_conf)
print('* profile length and azimuth range and interval: xlenp,alphab,alphae,dalp', file=binr_conf)
print("%.1f,%.1f,%.1f,%.1f" %
      (Profile_len, az_min, az_max, az_step), file=binr_conf)
print('* the spacing between bins, least number of traces, rnumtra, UTM_PROJECTION_ZONE(new)', file=binr_conf)
print("%d %d %.1f %d" % (bins_step, trace_num_min,
      ratio_trace, UTM_zone), file=binr_conf)
print('* time file name: timefile', file=binr_conf)
############
print("../m660q/%s" % m660q_Sp_out, file=binr_conf)
############
print('* output file name: outfile', file=binr_conf)
# print('inw20cw_ispwnccaz90_yb15-100vnt2_xb200_dx2_norm0_nf2p5-s1_Pcs',file=binr_conf)
print("%s" % binr_out_name, file=binr_conf)
print('* ouput number of time samples in each trace: npt, dt', file=binr_conf)
print("%d     %.2f" % (out_trace_npts, out_trace_dt), file=binr_conf)
print('* the indexes of reference ray among 1 -- nw: ninw, (inw0(i),i=1,ninw)   ', file=binr_conf)
print('20,1 -400,2 120,3 120,4 120,5 120,6 120,7 120,8 120,9 120,10 120,11 120,12 120,13 120,14 120,15 120,16 120,17 120,18 120,19 120,20 120', file=binr_conf)
# 1 -200 represent all depth with XBIN=200km
print('* minimum YBIN (km)', file=binr_conf)
print('15,20,25,30,30,32,32,34,34,36,36,38,38,40,40,42,44,46,48,50', file=binr_conf)
print('* DYBIN (km)', file=binr_conf)
print('-50,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,4', file=binr_conf)
# First -2 represent all depth with DYBIN=2km
print('* maximum YBIN (km)', file=binr_conf)
print('30,40,50,60,60,65,65,70,70,75,75,80,80,80,80,85,90,90,95,100', file=binr_conf)
print('* temporary directory name to store the intermedial files (.img)', file=binr_conf)
print("%s" % temp_folder, file=binr_conf)
print('* moveout index: idist, gcarc1  (only useful for idist=1)', file=binr_conf)
print("%d  %.2f" % (moveout_flag, moveout_gcarc), file=binr_conf)
print('* inorm', file=binr_conf)
print("%d " % norm_flag, file=binr_conf)
print('* output number and depth indexes in ninw: noutd,(ioutd(i),i=1,noutd)', file=binr_conf)
print('5 6 13 16 19 36      40-, 100-, 150-, 210-, 300-km', file=binr_conf)
print('* output index for stacking: istack', file=binr_conf)
print("%d " % stack_flag, file=binr_conf)
print('* output index for gcarc, baz and p: ioutb', file=binr_conf)
print("%d" % out_idx_flag, file=binr_conf)
print('* piercing point data file number: npief', file=binr_conf)
print("%d" % npief, file=binr_conf)
print('* input file name: infile', file=binr_conf)
print("%s" % pierc_out, file=binr_conf)
binr_conf.close()

# STEP 03
os.chdir(PSDM_pro_path+'stack')
if os.path.exists("./bin.out"):
    os.remove("./bin.out")
if os.path.exists("./bin.dat"):
    os.remove("./bin.dat")
if os.path.exists("./"+binr_out_name+"_num.dat"):
    os.remove("./"+binr_out_name+"_num.dat")
if os.path.exists("./"+binr_out_name+"_profile.txt"):
    os.remove("./"+binr_out_name+"_profile.txt")
if os.path.exists("./"+binr_out_name+"_yb.dat"):
    os.remove("./"+binr_out_name+"_yb.dat")
if os.path.exists("./stack_"+binr_out_name+".dat"):
    os.remove("./stack_"+binr_out_name+".dat")
print(">>> Step03 # Processing in folder:\n  %s\n\nbinr_vary_scan_n INFO :\n" % os.getcwd())
proc03_make = Popen(
    'make clean\nmake\n',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
proc03_make.communicate()
proc03 = Popen(
    './binr_vary_scan_n',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
outinfo03, errinfo03 = proc03.communicate()
# print(outinfo03.decode())
# print(errinfo03.decode())

##################################################################################################
hdpmig_conf = open((PSDM_pro_path+'poststack/hdpmig.in'), 'w+')
print('* imethod (phshift=0; phscreen=1, hybscreen: else),irefvel,vscale ', file=hdpmig_conf)
print("%d %d %.2f" % (method_flag, refvel_flag, vscale), file=hdpmig_conf)
print('* fmin, fmax (Minimum and maximum frequencies), ifreqindl, ifreqindr', file=hdpmig_conf)
print("%f %f %d %d" % (freq_min, freq_max,
      freq_idx_b, freq_idx_e), file=hdpmig_conf)
print('* nxmod, nzmod, nx, nz', file=hdpmig_conf)
print("%d %d %d %d" % (nxmod, nzmod, nx, nz), file=hdpmig_conf)
print('* dx, dz', file=hdpmig_conf)
print("%d %.1f" % (dx, dz), file=hdpmig_conf)
print('* ntrace, nt, dt (in sec.), nt0, ntb', file=hdpmig_conf)
print("%d %d %.1f %d %d" % (ntrace, nt, dt, nt0, nt_step), file=hdpmig_conf)
print('* FD method (15, 45, 65)', file=hdpmig_conf)
print("%d" % FD_ang, file=hdpmig_conf)
print('* nxleft, nxright', file=hdpmig_conf)
print("%d %d" % (smoothpt_left, smoothpt_right), file=hdpmig_conf)
print('* ifmat (=0: ascii vel. file; else: binary vel. file)', file=hdpmig_conf)
print("%d" % mod_format_flag, file=hdpmig_conf)
print('* modvelocity', file=hdpmig_conf)
print("%s" % mig_model, file=hdpmig_conf)
print('* tx_data (input seismic data)', file=hdpmig_conf)
print("%s" % (stacked_field), file=hdpmig_conf)
print('* migdata (output imaging data)', file=hdpmig_conf)
print("%s" % (hdpmig_out_field), file=hdpmig_conf)
# parameter use in debug
print('* intrace', file=hdpmig_conf)
print('1', file=hdpmig_conf)
print('* first trace index: itrfirst', file=hdpmig_conf)
print('1', file=hdpmig_conf)
hdpmig_conf.close()
# STEP 04
os.chdir(PSDM_pro_path+'poststack')
if os.path.exists("./"+hdpmig_out_field):
    os.remove("./"+hdpmig_out_field)
proc04_make = Popen(
    'make clean\nmake\n',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
proc04_make.communicate()
print(">>> Step04 # Processing in folder:\n  %s\n\nhdpmig.x INFO :\n" % os.getcwd())
proc04 = Popen(
    './hdpmig.x',
    stdin=None,
    stdout=PIPE,
    stderr=PIPE,
    shell=True)
outinfo04, errinfo04 = proc04.communicate()
# print(outinfo04.decode())
# print(errinfo04.decode())

print(errinfo02.decode())
print(errinfo03.decode())
print(errinfo04.decode())

# Observation System PLOT PART
##################################################################################################
config = {"font.family": 'sans-serif',
          "font.size": 16, "mathtext.fontset": 'stix'}
rcParams.update(config)
sta_xy_file = PSDM_pro_path+"/stack/station.dat"
pierce_xy_file1 = PSDM_pro_path+"/stack/depth1.dat"
pierce_xy_file2 = PSDM_pro_path+"/stack/depth2.dat"

pierce_xy_file3 = PSDM_pro_path+"/stack/depth3.dat"
pierce_xy_file4 = PSDM_pro_path+"/stack/depth4.dat"
pierce_xy_file5 = PSDM_pro_path+"/stack/depth5.dat"
pro_xy_file = PSDM_pro_path+"/stack/"+binr_out_name+"_profile.txt"
sta_xy = pd.read_csv(sta_xy_file, sep='\s+', header=None,
                     names=["lat", "lon", "decy", "decx"])
pierce_xy1 = pd.read_csv(pierce_xy_file1, sep='\s+',
                         header=None, names=["lat", "lon", "decy", "decx"])
pierce_xy2 = pd.read_csv(pierce_xy_file2, sep='\s+',
                         header=None, names=["lat", "lon", "decy", "decx"])
pierce_xy3 = pd.read_csv(pierce_xy_file3, sep='\s+',
                         header=None, names=["lat", "lon", "decy", "decx"])
pierce_xy4 = pd.read_csv(pierce_xy_file4, sep='\s+',
                         header=None, names=["lat", "lon", "decy", "decx"])
pierce_xy5 = pd.read_csv(pierce_xy_file5, sep='\s+',
                         header=None, names=["lat", "lon", "decy", "decx"])
pro_xy = pd.read_csv(pro_xy_file, sep='\s+', header=None,
                     names=["lon", "lat", "decx", "decy", "offset"])
# CD_xy=pd.read_csv("/home/georom1996/PracticeJoe/PyPSDM/SW_China_Vs_model/sw.lst",sep='\s+',header=None,names=["lon","lat"])
plt.figure(figsize=(20, 8))
# plt.set_title("Observation System View")
# plt.scatter(pierce_xy.lon,pierce_xy.lat)
# plt.plot(pro_xy.lon,pro_xy.lat)
plt.subplot(121)
# plt.scatter(pierce_xy3.decx,pierce_xy3.decy,color='gray',marker='.')
plt.scatter(pierce_xy2.decx, pierce_xy2.decy, color='green', marker='.')
plt.scatter(pierce_xy1.decx, pierce_xy1.decy, color='red', marker='.')

# plt.scatter(pierce_xy4.decx,pierce_xy4.decy,color='brown',marker='+')
# plt.scatter(pierce_xy5.decx,pierce_xy5.decy,color='cyan',marker='+')
plt.scatter(sta_xy.decx, sta_xy.decy, color='k', marker='^',)
# plt.scatter(CD_xy.lon,CD_xy.lat,color='blue',marker='^')
plt.plot(pro_xy.decx, pro_xy.decy, 'blue', linewidth=3)

# plt.axis('equal')
plt.subplot(122)

# plt.scatter(pierce_xy3.lon,pierce_xy3.lat,color='gray',marker='.')
plt.scatter(pierce_xy2.lon, pierce_xy2.lat, color='green', marker='.')
plt.scatter(pierce_xy1.lon, pierce_xy1.lat, color='red', marker='.')
# plt.scatter(pierce_xy4.lon,pierce_xy4.lat,color='brown',marker='+')
# plt.scatter(pierce_xy5.lon,pierce_xy5.lat,color='cyan',marker='+')
plt.scatter(sta_xy.lon, sta_xy.lat, color='k', marker='^')
# plt.scatter(CD_xy.lon,CD_xy.lat,color='blue',marker='^')
plt.plot(pro_xy.lon, pro_xy.lat, 'blue', linewidth=3)
plt.savefig(cwdfd+"/figs/"+"example002_IASP_Station_Map.jpg",
            dpi=300, bbox_inches='tight')

# plt.axis('equal')

# CCP PLOT PART
##################################################################################################
# Plot CCP stacking section with bin size and RF numbers at chosen depth
# Inputs (parameters set in binr_vary_scan_n.inp):
# 1. outfile    -   output file name, no prefix and postfix
#                   e.g. stack_*.dat, *_yb.dat, *_num.dat
# 2. xlenp      -   profile length (km)
# 3. npt        -   output number of time samples in each trace
# 4. dx         -   the spacing between bins along the profile (km)
# 5. dt         -   time interval (s)
# 6. noutd      -   output number in ninw
config = {"font.family": 'sans-serif',
          "font.size": 22, "mathtext.fontset": 'stix'}
rcParams.update(config)


def add_right_cax(ax, pad, width):
    '''
    在一个ax右边追加与之等高的cax.
    pad是cax与ax的间距.
    width是cax的宽度.    '''
    axpos = ax.get_position()
    caxpos = mpl.transforms.Bbox.from_extents(
        axpos.x1 + pad,
        axpos.y0,
        axpos.x1 + pad + width,
        axpos.y1
    )
    cax = ax.figure.add_axes(caxpos)
    return cax


##################################################################################################
xlenp = 1000
npt = 1001
dt = 0.1
depth = 800  # imaging depth (km)
dz = 0.5  # depth interval, km per grid along Z
prolen = 1000  # profile length
dx = 50  # trace interval, km per grid along X
nz = int(depth/dz)  # depth number
nx = int(prolen/dx)+1  # trace number
Yrange = np.linspace(0, nz*dz, nz)
Trange = np.linspace(0, npt*dt, npt)
Xrange = np.linspace(0, xlenp, int(xlenp/dx)+1)
noutd = 5
numbin = int(xlenp/dx+1)
timelen = (npt-1)*dt
stacked_field = PSDM_pro_path+"/stack/stack_"+binr_out_name+".dat"
yb_file = PSDM_pro_path+"/stack/"+binr_out_name+"_yb.dat"
num_file = PSDM_pro_path+"/stack/"+binr_out_name+"_num.dat"
LatRange = pro_xy.lat
with open(stacked_field, 'rb') as f:
    for k in range(1):
        data = np.fromfile(f, dtype=np.float32)
        profile = (np.reshape(data, (numbin, npt)))
with open(yb_file, 'rb') as f_yb:
    for k_yb in range(1):
        data_yb = np.fromfile(f_yb, dtype=np.float32)
        yb_data = (np.reshape(data_yb, (noutd, numbin)))
with open(num_file, 'rb') as f_num:
    for k_yb in range(1):
        data_num = np.fromfile(f_num, dtype=np.float32)
        num_data = (np.reshape(data_num, (noutd, numbin)))

method = "Qseis"
vol_data = "Mig"
filename = PSDM_pro_path+"/poststack/"+hdpmig_out_field
with open(filename, 'rb') as f:
    for k in range(1):
        data_mig = np.fromfile(f, dtype=np.float32)
        profile_mig = (np.reshape(data_mig, (nx, nz)))

fig, ax = plt.subplots(3, 1, figsize=(16, 14), sharex=True,
                       gridspec_kw={'height_ratios': [2, 3, 3]})
im00 = ax[0].plot(LatRange, yb_data[1], 'k', lw=4, label="Half Bin Width(km)")
ax[0].set_ylabel("Half Bin Width(km)")
# ax[0].legend(["Half Bin Width (km)"])
ax_00 = ax[0].twinx()
ax_00.plot(LatRange, num_data[1, :], 'r', label="RF num")
# ax_00.set_ylim([0,400])
ax_00.set_ylabel("RF number")
ax_00.set_title("Parameters at 100km")
ax[0].legend(loc=2)
ax_00.legend(loc=1)

im1 = ax[1].pcolor(LatRange, Trange, profile.T,
                   cmap="coolwarm", vmin=-0.6, vmax=0.6)
# im1=ax[1].contourf(Xrange,Trange,profile.T)
cax = add_right_cax(ax[1], pad=0.02, width=0.02)
cbar = fig.colorbar(im1, cax=cax)
im1.set_clim(-0.5, 0.5)
ax[1].set_ylim([0, 30])
# ax[1].set_xlabel("Offset (km)")
ax[1].set_ylabel("T(s)")
ax[1].invert_yaxis()
ax[1].set_title("CCP stack section")
# fig.colorbar()
im1 = ax[2].pcolormesh(LatRange, Yrange, profile_mig.T,
                       cmap="jet", vmin=-0.6, vmax=0.6)
# im1=ax[1].contourf(Xrange,Trange,profile.T)
cax = add_right_cax(ax[2], pad=0.02, width=0.02)
cbar = fig.colorbar(im1, cax=cax)
im1.set_clim(-0.45, 0.55)
ax[2].set_ylim([0, 250])
ax[2].set_xlabel("Latitute ($\degree$)")
ax[2].set_ylabel("Depth(km)")
ax[2].grid()
ax[2].plot(stalstdata.stla, stalstdata.moho,
           color='white', marker='o', linestyle='', ms=5)
ax[2].plot(stalstdata.stla, stalstdata.lab,
           color='white', marker='o', linestyle='', ms=5)
ax[2].plot(stalstdata.stla, 5+np.zeros(len(stalstdata.stla)),
           color='r', linestyle='', marker='v', ms=10)
ax[2].invert_yaxis()
ax[2].set_title(method+"_"+str(vol_data))
ax[2].text(0.65, 0.1, (method+"_"+str(vol_data)+" Freq "+str(freq_min)+"-" +
           str(freq_max)+"Hz"), fontdict={'size': '24', 'color': 'k'}, transform=ax[2].transAxes)
# ax[0].set_title("Migration section")
# ax[2].text(0.65,0.02,(method+"_"+str(vol_data)+" Freq "+str(freq_min)+"-"+str(freq_max)+"Hz"),fontdict={'size':'24','color':'k'},transform=ax[2].transAxes)
# plt.tight_layout()
# plt.savefig("/home/georom1996/PracticeJoe/PyPSDM/SynPSDM/QseisMigration2022.png",bbox_inches='tight')
# plt.savefig("/home/georom1996/PracticeJoe/PyPSDM/SynPSDM/"+model_tag+"QseisMigration2022.pdf",format='pdf')
plt.savefig(cwdfd+"/figs/"+"example002_IASP_PSDM_Map.jpg",
            dpi=300, bbox_inches='tight')
