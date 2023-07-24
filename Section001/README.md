## 操作流程
- 进入src文件夹,编译相关源码
```bash
cd src
make clean
make all
cd ../
```   

- 进入example01_kumar文件夹，执行复现Kumar文章的案例
```bash
cd ./example01_kumar
bash ./SX_Deconv_example_kumar.sh
cd ../
```
结果放置于example01_kumar/figs文件夹中，图件如下：   
![Example Kumar](./example01_kumar/figs/corrected-XR.ST09.1998.197.121411-SRF.l.jpg)

- 进入example02_gcsrf文件夹，执行GC_SRF计算S波接收函数的案例
```bash
cd ./example02_gcsrf
bash ./GC_SRF_example_syn.sh
```   
然后执行绘图命令
```bash
python ./GC_SRF_plot.py
```   
结果图件放置于example02_gcsrf/figs文件夹中，图件如下：   
![GC_SRF 001](./example02_gcsrf/figs/syn_SM_d75_001.jpg)   
![GC_SRF 002](./example02_gcsrf/figs/syn_SM_d75_002.jpg)

## 必要工具和库
### Linux
1. Taup
2. GMT 6
3. SAC
4. MPIC++
5. GCC

### Python Lib
1. obspy 
2. pandas 
3. numpy 
4. matplotlib