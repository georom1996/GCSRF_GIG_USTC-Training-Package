# Function Part : JoeZhang 2021/08/12 Email:zhangzhou1996@gmail.com
import os
from obspy import read
import pandas as pd
# from pandas import read_csv,read_csv
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

def autoscale_x(ax, margin=0.1):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

    import numpy as np

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo, hi = ax.get_ylim()
        x_displayed = xd[((yd > lo) & (yd < hi))]
        h = np.max(x_displayed) - np.min(x_displayed)
        bot = np.min(x_displayed)-margin*h
        top = np.max(x_displayed)+margin*h
        return bot, top

    lines = ax.get_lines()
    bot, top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot:
            bot = new_bot
        if new_top > top:
            top = new_top

    ax.set_xlim(bot, top)


def autoscale_y(ax, margin=0.1):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

    import numpy as np

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo, hi = ax.get_xlim()
        y_displayed = yd[((xd > lo) & (xd < hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot, top

    lines = ax.get_lines()
    bot, top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot:
            bot = new_bot
        if new_top > top:
            top = new_top

    ax.set_ylim(bot, top)


def GCSRF_plotSyn(SRF_fd, fd_lst, fd_begin, ta, tb, tc, td, te, oname):
    for i in fd_lst:
        if (i.startswith(fd_begin)):
            # if(i==fdname):
            print(i)
            fdname = i
            # amp=pd.read_csv(SRF_fd+"/"+fdname+"/all.Amp_0.xyz",sep=' ',header=None, names=['inci','winlen','value'])

            x, y, z = np.loadtxt(
                SRF_fd+"/"+fdname+"/all.Amp_0.xyz", unpack=True)
            nx = len(np.unique(x))
            ny = len(np.unique(y))
            x1, y1, z1 = np.loadtxt(
                SRF_fd+"/"+fdname+"/all.CCC.xyz", unpack=True)
            x2, y2, z2 = np.loadtxt(
                SRF_fd+"/"+fdname+"/all.Energy_0.xyz", unpack=True)
            # min1,min2,min3,min4=np.loadtxt(SRF_fd+"/"+fdname+"/Min.lst",unpack=True)
            Minlst = pd.read_csv(SRF_fd+"/"+fdname+"/Min.lst", sep=' ',
                                header=None, names=['inci', 'winlen', 'amp', 'name'])
            Synlst = pd.read_csv(SRF_fd+"/"+fdname+"/Syn.dat",
                                sep=' ', header=None, names=['col1', 'col2', 'inci'])
            OptSac = pd.read_csv(SRF_fd+"/"+fdname+"/OptimalSRF.lst", sep=' ',
                                header=None, names=['name', 'inci', 'win_len', 'cc'])

            plt.figure(figsize=(5, 8))
            # plt.rcParams['figure.figsize']=(5.0,12.0)
            plt.rcParams['figure.dpi'] = 100
            plt.subplot(2, 1, 1)
            # plt.contourf(x.reshape(ny,nx),y.reshape(ny,nx),z.reshape(ny,nx),50,cmap=plt.cm.bwr,vmin=-1,vmax=1)

            plt.pcolormesh(x.reshape(ny, nx), y.reshape(ny, nx), z.reshape(
                ny, nx), cmap=plt.cm.coolwarm, vmin=-1, vmax=1)
            plt.colorbar()

            plt.plot(Minlst.inci, Minlst.winlen, '+', color='g')
            plt.plot(float(OptSac.name[0].split(".")[1]), float(
                OptSac.name[0].split(".")[2]), '+', color='r', lw=2)
            plt.plot([Synlst.inci]*len(Minlst.winlen),
                     Minlst.winlen, '--', color='k')
            plt.xlabel('Inci_Ang (degree)')
            plt.ylabel('Win_Len (s)')
            plt.title('Amp_0')
            plt.text(2, 10, ta, fontsize=18, bbox=dict(
                facecolor='white', alpha=0.8))
            # plt.legend(["(a)"])

            plt.subplot(2, 1, 2)
            plt.pcolormesh(x1.reshape(ny, nx), y1.reshape(ny, nx), z1.reshape(
                ny, nx), cmap=plt.cm.plasma, vmin=0, vmax=1)
            plt.plot(Minlst.inci, Minlst.winlen, '+', color='g')
            plt.plot(float(OptSac.name[0].split(".")[1]), float(
                OptSac.name[0].split(".")[2]), '+', color='r', lw=2)
            plt.plot([Synlst.inci]*len(Minlst.winlen),
                     Minlst.winlen, '--', color='k')
            plt.xlabel('Inci_Ang (degree)')
            plt.ylabel('Win_Len (s)')
            plt.title('Coef_SRF')
            plt.text(2, 10, tb, fontsize=18, bbox=dict(
                facecolor='white', alpha=0.8))
            # plt.legend(["(b)"])
            plt.colorbar()
            plt.tight_layout()
            plt.savefig("./figs/"+fdname+'_001.pdf', dpi=600, format='PDF')
            plt.savefig("./figs/"+fdname+'_001.jpg')
            # plt.figure(figsize=(6,6))
            fig = plt.figure(figsize=(8, 8), dpi=100)

            grid = plt.GridSpec(4, 4, hspace=0.5, wspace=0.5)
            main_ax = fig.add_subplot(grid[1:4, 0:3])
            rcc_ax = fig.add_subplot(grid[1:4, 3:4], sharey=main_ax)
            tss_ax = fig.add_subplot(grid[0:1, 0:3], sharex=main_ax)
            # plt.subplot(1,4,4)
            RedSac = pd.read_csv(SRF_fd+"/"+fdname+"/SelectedRed.CC.lst",
                                sep=' ', header=None, names=['name', 'inci', 'win_len', 'cc'])
            OptSac = pd.read_csv(SRF_fd+"/"+fdname+"/OptimalSRF.lst", sep=' ',
                                header=None, names=['name', 'inci', 'win_len', 'cc'])
            mul = 20
            allsrf = sorted(os.listdir(SRF_fd+"/"+fdname+"/"))
            sacend = ".cutrv"
            for srf in allsrf:
                if (srf.endswith(sacend)):
                    st = read(SRF_fd+"/"+fdname+"/"+srf)
                    st.normalize()
                    win_len = int(srf.split('.')[2])
                    tr = st[0]
                    begin = st[0].stats.sac.b
                    end = st[0].stats.sac.e
                    dt = st[0].stats.sac.delta
                    data = st[0].data
                    time = np.linspace(begin, end, len(data))
                    main_ax.plot(time, data*mul+win_len, color='lightgray')
                    main_ax.set_xlabel('Time (s)')
                    main_ax.set_ylabel('Win_Len (s)')
                    main_ax.set_xlim([0, 30])
                    main_ax.set_ylim([0, 100])
            main_ax.text(0.75, 4, td, fontsize=18,
                         bbox=dict(facecolor='white', alpha=0.8))

            sumdata = np.zeros(len(data))
            idxnum = 0
            for idx in RedSac.index:
                st = read(SRF_fd+"/"+fdname+"/" +
                          str(RedSac.name[idx].replace('cutrv.rm.rv', 'cutrv', 1)))
                st.normalize()
                win_len = RedSac.win_len[idx]
                tr = st[0]
                begin = st[0].stats.sac.b
                end = st[0].stats.sac.e
                dt = st[0].stats.sac.delta
                data = st[0].data
                sumdata += data
                time = np.linspace(begin, end, len(data))
                main_ax.plot(time, data*mul+win_len, color='g')
                main_ax.set_xlim([0, 30])
                main_ax.set_ylim([0, 100])
                idxnum += 1
            tss_lg = []

            # DoRF_Path="/home/georom1996/PracticeJoe/RF_modeling.joe/DoSynRF.joe/"
            # test_tag="TestMod_20211022"
            # if oname == "SM01":
            #     modname="MOD03"
            # elif oname == "SM02":
            #     modname="MOD01"
            # elif oname == "SM03":
            #     modname="MOD02"

            # dosyn_path=DoRF_Path+test_tag+"/"+modname+"-srf.dd"
            # srfdd11=pd.pd.read_csv(dosyn_path,header=None,names=["data"])
            # srfdd=np.flip(srfdd11.data.to_numpy())
            # plt.plot(srfdd)
            nfft = 2048
            dt = 0.1
            time_srf = np.linspace(100-nfft*dt, 100, nfft)
            # tss_ax.plot(time_srf,-8*srfdd,'c',lw=4)
            # tss_lg.append("SynSRF")
            tss_ax.text(0.75, 0.75, tc, fontsize=18,
                        bbox=dict(facecolor='white', alpha=0.8))
            # autoscale_y(tss_ax)

            rcc_ax.plot(RedSac.cc, RedSac.win_len, 'o', color='g')
            for idx in OptSac.index:
                st = read(SRF_fd+"/"+fdname+"/" +
                          str(OptSac.name[idx].replace('cutrv.rm.rv', 'cutrv', 1)))
                st.normalize()
                win_len = int(OptSac.name[idx].split('.')[2])
                tr = st[0]
                begin = st[0].stats.sac.b
                end = st[0].stats.sac.e
                dt = st[0].stats.sac.delta
                data = st[0].data
                time = np.linspace(begin, end, len(data))
                main_ax.plot(time, data*mul+win_len, color='r')
                tss_ax.plot(time, data, 'r', lw=4)
                tss_lg.append("OptSRF")

                main_ax.set_xlim([0, 30])
                main_ax.set_ylim([0, 105])
                rcc_ax.plot(OptSac.cc, win_len, 'o', color='r')

                # rcc_ax.arrow(0.1, win_len, OptSac.cc-0.1, 0.5,length_includes_head=True,head_width=0.25, head_length=0.5, fc='r', ec='b')
                # rcc_ax.quiver(OptSac.cc-0.1,,win_len,10, 0, color='g', width=0.1)
                rcc_ax.annotate('OptSRF', xy=(OptSac.cc, win_len), xytext=(
                    OptSac.cc-0.8, win_len-3), color='k', arrowprops=dict(arrowstyle='->'))
            # print(OptSac)
            tss_ax.plot(time, sumdata/idxnum, 'k--', lw=2)
            rcc_ax.text(0.25, 98, te, fontsize=18,
                        bbox=dict(facecolor='white', alpha=0.8))
            rcc_ax.set_title("Coef_SRF")
            tss_lg.append("SumAll")
            # print(tss_lg)
            tss_ax.legend(tss_lg)
            fig.savefig("./figs/"+fdname+'_002.pdf', dpi=600, format='PDF')
            fig.savefig("./figs/"+fdname+'_002.jpg')
#                 fig.tight_layout()
#                 fig.show()

            # plt.figure()
            # plt.plot(srfdd)

            # plt.pcolor(x.reshape(ny,nx),y.reshape(ny,nx),z.reshape(ny,nx),50,cmap=plt.cm.seismic)
            # ccc=pd.read_csv(SRF_fd+"/"+fdname+"/all.CCC.xyz",sep=' ',header=None, names=['inci','winlen','value'])
            # ene=pd.read_csv(SRF_fd+"/"+fdname+"/all.Energy_0.xyz",sep=' ',header=None, names=['inci','winlen','value'])
