import pandas as pd
# import obspy
from obspy.taup.taup_create import build_taup_model
mod_text=pd.read_csv('./CDmod',sep="\s+",header=0,names=['vp','vs','rho','thick','depth','idx'])
iasp_text=pd.read_csv('./iasp91.tvel',sep="\s+",header=2,names=['depth','vp','vs','rho'])
fileout=open("./CDmod.tvel","w")
print( "CDmod P Model (%d layers) no \"discontinuity\" at 120, 760 km" % (len(mod_text)),file=fileout)
print( "CDmod S Model (%d values to cmb)" % (len(mod_text)+iasp_text[iasp_text.depth==2939.330].index.tolist()[0]),file=fileout)
for ii in range(len(mod_text)):
    print(mod_text.depth[ii],mod_text.vp[ii],mod_text.vs[ii],mod_text.rho[ii],file=fileout)
idx_0=iasp_text[iasp_text.depth>mod_text.depth[ii]].index.tolist()[0]
for jj in range(idx_0,len(iasp_text)):
    print(iasp_text.depth[jj],iasp_text.vp[jj],iasp_text.vs[jj],iasp_text.rho[jj],file=fileout)
fileout.close()
build_taup_model("./CDmod.tvel")