## 操作流程
- 进入src文件夹,编译相关源码
```bash
cd src
make clean
make all
cd ../
```   
- 进入example01_kumar文件夹,执行复现Kumar文章的案例
```bash
cd ./example01_kumar
bash ./SX_Deconv_example_kumar.sh
cd ../
```
结果放置于example01_kumar/figs文件夹中,图件如下：   
![Example Kumar](./example01_kumar/figs/corrected-XR.ST09.1998.197.121411-SRF.l.jpg)

- 进入example02_gcsrf文件夹,执行GC_SRF计算S波接收函数的案例
```bash
cd ./example02_gcsrf
bash ./GC_SRF_example_syn.sh
```   
然后执行绘图命令
```bash
python ./GC_SRF_plot.py
```   
结果图件放置于example02_gcsrf/figs文件夹中,图件如下：   
![GC_SRF 001](./example02_gcsrf/figs/syn_SM_d75_001.jpg)   
![GC_SRF 002](./example02_gcsrf/figs/syn_SM_d75_002.jpg)

## 文件说明
- 在流程操作编译完成之后,在bin文件夹中将出现如下命令：
    - 1. **SX_Deconv**
            该程序为脉冲反褶积计算的独立计算程序,使用方法为：
            ```bash
            SX_Deconv -Z Z.sac -N N.sac -E E.sac –reg 3 –inc 16 -t1 100 -t2 120 -D 
            ```
            程序选项：
            ```
                -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间)
                -t2 脉冲结束时刻 (相对SAC文件开始时刻的时间)
                -reg 计算逆滤波器时对角元素系数
                -inc 旋转角度
                -D 输出逆滤波器及波形旋转后的SAC文件
            ```   
    - 2. **GC_SRF**
            该程序为循环不同子波时窗长度和旋转角度循环计算反褶积的独立程序,
            程序使用方法为：
            ```bash
            GC_SRF -Z Z.sac -N N.sac -E E.sac -t1 t1 –reg 3 
            ```
            程序选项：
            ```
                -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间)
                -reg 计算逆滤波器时对角元素系数
                -D 输出逆滤波器及波形旋转后的SAC文件  

            程序中默认循环 win_len=0:5:100; inci=0:4:60
            其中t1一般设定为直达S到时相对SAC文件起始时刻的时间, win_len起算点为t1-5s
            
            ```
    - 3. **GC_SRF_par**
            该程序为GC_SRF的并行版本,使用方法与GC_SRF相同,用于加快循环计算的效率：
            ```bash
            GC_SRF_par -Z Z.sac -N N.sac -E E.sac -t1 t1 –reg 3 
            ```
            程序选项：
            ```
                -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间)
                -reg 计算逆滤波器时对角元素系数
                -D 输出逆滤波器及波形旋转后的SAC文件
            ```
    - 4. **GC_SRF_post**
            该程序为GC_SRF/GC_SRF_par计算完成的后处理程序,用于计算所有反褶积结果波形的0时刻振幅能量网格,互相关系数网格,并：
            ```bash
            GC_SRF_par -Z Z.sac -N N.sac -E E.sac -t1 t1 –reg 3 
            ```
            程序选项：
            ```
                -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间)
                -reg 计算逆滤波器时对角元素系数
                -D 输出逆滤波器及波形旋转后的SAC文件
            ```
    - 5. spmove.joe
    - 6. sacsnr
    - 7. sacStack
- 其中1~4为计算S波接收函数的主要程序,spmove.joe修改自Xiaohui Yuan教授spmout代码。sacanr和sacStack来自[SACTOOLS](https://github.com/msthorne/SACTOOLS/tree/master).


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

