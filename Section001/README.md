# README

本仓库提供中英文两种语言的文档说明。

- [English](#english-version)
- [Chinese](#chinese-version)

## English Version <a name="english-version"></a>

### Section001 Implementation of Pulse Deconvolution and GC_SRF Strategy

#### Basic Program Compilation

- Enter the `src` folder and compile the relevant source code:

```bash
cd src
make clean
make all
cd ../
```

#### Explanation of Executable Programs

After compilation, the following executable programs will appear in the `bin` folder:

1. **SX_Deconv**

   This program is an independent calculation program for pulse deconvolution. Usage:

   ```bash
   SX_Deconv -Z Z.sac -N N.sac -E E.sac –inc 16 -t1 100 -t2 120 -D 
   ```

   Program options:

   ```
   -t1 Pulse start time (relative to the start time of the SAC file)
   -t2 Pulse end time (relative to the start time of the SAC file)
   -inc Rotation angle
   -D Output inverse filter and waveform after rotation SAC file
   ```

   Input files:

   ```
   Z.sac N.sac E.sac
   Note: Ensure that baz exists in the header
   ```

   Output files:

   ```
   Default output:
   rf.sl.sac : Unnormalized L-component deconvolved waveform (SRF)
   rf.sq.sac : Unnormalized Q-component deconvolved waveform
   rf.st.sac : Unnormalized T-component deconvolved waveform
   rf.nsl.sac : Normalized L-component deconvolved waveform (SRF)
   rf.nsq.sac : Normalized Q-component deconvolved waveform
   rf.nst.sac : Normalized T-component deconvolved waveform
   If -D option is used, additionally output:
   l-com.sac : Original L-component waveform
   q-com.sac : Original Q-component waveform
   t-com.sac : Original T-component waveform
   fil.sac : Inverse filter
   ```

2. **GC_SRF**

   This program is an independent program for looping different subwave time window lengths and rotation angles to calculate deconvolutions. Usage:

   ```bash
   GC_SRF -Z Z.sac -N N.sac -E E.sac -t1 t1
   ```

   Program options:

   ```
   -t1 Pulse start time (relative to the start time of the SAC file)
   Note: By default, the program loops win_len=0:5:100; inci=0:4:60, where t1 is set as the arrival time of the direct S wave relative to the start time of the SAC file, and win_len starts from t1-5s
   ```

   Input files:

   ```
   Z.sac N.sac E.sac
   Note: Ensure that baz exists in the header
   ```

   Output files:

   ```
   Default output:
   ("rf-%02.f-%03.f.sl.sac", incnum, win_len) : Unnormalized L-component deconvolved waveform (SRF)
   ("rf-%02.f-%03.f.sq.sac", incnum, win_len) : Unnormalized Q-component deconvolved waveform
   ("rf-%02.f-%03.f.st.sac", incnum, win_len) : Unnormalized T-component deconvolved waveform
   ("rf-%02.f-%03.f.nsl.sac", incnum, win_len) : Normalized L-component deconvolved waveform (SRF)
   ("rf-%02.f-%03.f.nsq.sac", incnum, win_len) : Normalized Q-component deconvolved waveform
   ("rf-%02.f-%03.f.nst.sac", incnum, win_len) : Normalized T-component deconvolved waveform
   ```

3. **GC_SRF_par**

   This is the parallel version of GC_SRF, with the same usage as GC_SRF, used to speed up the looping calculations:

   ```bash
   GC_SRF_par -Z Z.sac -N N.sac -E E.sac -t1 t1
   ```

   Program options:

   ```
   -t1 Pulse start time (relative to the start time of the SAC file)
   -D Output inverse filter and waveform after rotation SAC file
   ```

4. **GC_SRF_post**

   This is the post-processing program for GC_SRF/GC_SRF_par, used to calculate the zero-time amplitude energy grid and cross-correlation coefficient grid of all deconvolved waveforms:

   ```bash
   GC_SRF_post
   ```

   Program options:

   ```
   No options, run directly in the GC_SRF calculation directory
   ```

   Input files:

   ```
   SAC format data files ending with .nsl.sac.100.cutrv.rm.rv after processing by GC_SRF/GC_SRF_par

   Note: The input SAC format data files are obtained as follows:
   1. Cut the output waveforms from GC_SRF/GC_SRF_par (files ending with nsl.sac) from -100s to 100s (.100.cut);
   2. Reverse the time axis (rv);
   3. Apply moveout correction (.rm);
   4. Reverse the time axis again (rv).
   The specific operation script can be referred to in Section001/example02_gcsrf/GC_SRF_example_syn.sh
   ```

   Output files:

   ```
   all.Amp_0.xyz : Zero-time amplitude spectrum of all SRFs
   all.Energy_0.xyz : Zero-time energy spectrum of all SRFs
   Min.lst : List of waveforms with the lowest zero-time amplitude for different time window lengths, format: rotation angle, time window length, filename
   all.CCC.xyz : Cross-correlation coefficient spectrum of all SRFs with the mean waveform
   SelectedRed.CC.lst : List of cross-correlation coefficients corresponding to the waveforms with the lowest zero-time amplitude for different time window lengths, format: rotation angle, time window length, filename
   OptimalSRF.lst : Optimal SRF, format: filename, rotation angle, time window length, cross-correlation coefficient
   ```

5. **spmove.joe**

   Modified from Professor Xiaohui Yuan's spmout code

6. **sacsnr** and **sacStack**

   From [SACTOOLS](https://github.com/msthorne/SACTOOLS/tree/master)

#### Case Test 01: Reproducing Kumar et al., 2006 Figure 4b

- Enter the `example01_kumar` folder and run the test script `SX_Deconv_example_kumar.sh`:

```bash
cd ./example01_kumar
bash ./SX_Deconv_example_kumar.sh
cd ../
```

The results are placed in the `example01_kumar/figs` folder, as shown in the figure:

![Example Kumar](./example01_kumar/figs/corrected-XR.ST09.1998.197.121411-SRF.l.jpg)

The original figure from Kumar et al., 2006 Figure 4 can be seen at [Kumar et al., 2006](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JB003930)

#### Case Test 02: GC_SRF Calculation Process Based on Pulse Deconvolution

- Enter the `example02_gcsrf` folder and run the test script `GC_SRF_example_syn.sh`:

This script uses GC_SRF and GC_SRF_post programs to provide a reference process for calculating S-wave receiver functions based on GC_SRF.

```bash
cd ./example02_gcsrf
bash ./GC_SRF_example_syn.sh
```

- Then, use Python to run the plotting script:

```bash
python ./GC_SRF_plot.py
```

The result figures are placed in the `example02_gcsrf/figs` folder, as shown in the figures:

![GC_SRF 001](./example02_gcsrf/figs/syn_SM_d75_001.jpg)

![GC_SRF 002](./example02_gcsrf/figs/syn_SM_d75_002.jpg)

#### Necessary Tools and Libraries

##### Linux

1. Taup
2. GMT 6
3. SAC
4. MPIC++
5. GCC

##### Python Lib

1. obspy
2. pandas
3. numpy
4. matplotlib

## Chinese Version <a name="chinese-version"></a>

### Section001 脉冲反褶积与GC_SRF策略的实现

#### 基础程序编译

- 进入`src`文件夹,编译相关源码

```bash
cd src
make clean
make all
cd ../
```

#### 可执行程序说明

在流程操作编译完成之后,在`bin`文件夹中将出现如下可执行程序：

- 1. **SX_Deconv**

        该程序为脉冲反褶积计算的独立计算程序,使用方法为：

        ```bash
        SX_Deconv -Z Z.sac -N N.sac -E E.sac –inc 16 -t1 100 -t2 120 -D 
        ```

        程序选项：

        ```
            -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间)
            -t2 脉冲结束时刻 (相对SAC文件开始时刻的时间)
            -inc 旋转角度
            -D 输出逆滤波器及波形旋转后的SAC文件
        ```   

        输入文件：

        ```
            Z.sac N.sac E.sac
            注:头端中需要确保baz存在
        ```

        输出文件：

        ```
        默认输出：
            rf.sl.sac : 未归一化的L分量反褶积波形(SRF)
            rf.sq.sac : 未归一化的Q分量反褶积波形
            rf.st.sac : 未归一化的T分量反褶积波形
            rf.nsl.sac : 振幅归一化的L分量反褶积波形(SRF)
            rf.nsq.sac : 振幅归一化的Q分量反褶积波形
            rf.nst.sac : 振幅归一化的T分量反褶积波形
        若使用-D选项，则会额外输出：
            l-com.sac : L分量原始波形
            q-com.sac : Q分量原始波形
            t-com.sac : T分量原始波形
            fil.sac : 反子波滤波器
        ```

- 2. **GC_SRF**

        该程序为循环不同子波时窗长度和旋转角度循环计算反褶积的独立程序,
        程序使用方法为：

        ```bash
        GC_SRF -Z Z.sac -N N.sac -E E.sac -t1 t1
        ```

        程序选项：

        ```
            -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间) 
            注:程序中默认循环 win_len=0:5:100; inci=0:4:60, 其中 t1 设定为直达 S 到时相对 SAC 文件起始时刻的时间, win_len 起算点为 t1-5s
        ```

        输入文件：

        ```
            Z.sac N.sac E.sac
            注:头端中需要确保baz存在
        ```

        输出文件：

        ```
        默认输出：
            ("rf-%02.f-%03.f.sl.sac", incnum, win_len) : 未归一化的L分量反褶积波形(SRF)
            ("rf-%02.f-%03.f.sq.sac", incnum, win_len) : 未归一化的Q分量反褶积波形
            ("rf-%02.f-%03.f.st.sac", incnum, win_len) : 未归一化的T分量反褶积波形
            ("rf-%02.f-%03.f.nsl.sac", incnum, win_len) : 振幅归一化的L分量反褶积波形(SRF)
            ("rf-%02.f-%03.f.nsq.sac", incnum, win_len) : 振幅归一化的Q分量反褶积波形
            ("rf-%02.f-%03.f.nst.sac", incnum, win_len) : 振幅归一化的T分量反褶积波形
        ```

- 3. **GC_SRF_par**

        该程序为GC_SRF的并行版本,使用方法与GC_SRF相同,用于加快循环计算的效率：

        ```bash
        GC_SRF_par -Z Z.sac -N N.sac -E E.sac -t1 t1
        ```

        程序选项：

        ```
            -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间)
            -D 输出逆滤波器及波形旋转后的SAC文件
        ```

- 4. **GC_SRF_post**

        该程序为GC_SRF/GC_SRF_par的后处理程序,用于计算所有反褶积结果波形的0时刻振幅能量网格,互相关系数网格：

        ```bash
        GC_SRF_post
        ```

        程序选项：

        ```
            无任何选项，直接在GC_SRF计算目录执行
        ```

        输入文件：

        ```
            GC_SRF/GC_SRF_par输出波形处理后的以nsl.sac.100.cutrv.rm.rv结尾的SAC格式数据文件

            注:该输入SAC格式数据文件获取方式为：
            1. 将GC_SRF/GC_SRF_par输出波形(nsl.sac结尾文件)截取-100s到100s(.100.cut);
            2. 翻转时间轴(rv);
            3. moveout校正(.rm);
            4. 翻转转时间轴(rv).
            具体操作脚本可参考Section001/example02_gcsrf/GC_SRF_example_syn.sh
        ```

        输出文件：

        ```
            all.Amp_0.xyz : 所有SRF的0时刻振幅谱
            all.Energy_0.xyz : 所有SRF的0时刻能量谱
            Min.lst : 不同时窗长度对应的0时刻振幅最低的波形列表，文件格式为：旋转角度,时窗长度,文件名
            all.CCC.xyz : 所有SRF的与均值波形的互相关系数谱
            SelectedRed.CC.lst : 不同时窗长度0时刻振幅最低的波形对应的互相关系数列表，文件格式为：旋转角度,时窗长度,文件名
            OptimalSRF.lst : 最佳SRF,文件格式为：文件名,旋转角度,时窗长度,互相关系数
        ```

- 5. **spmove.joe**

        ```
            修改自Xiaohui Yuan教授spmout代码
        ```

- 6. **sacsnr** 和 **sacStack**

    来自[SACTOOLS](https://github.com/msthorne/SACTOOLS/tree/master).

#### 案例测试01： 复现kumar et al., 2006 Figure 4b

- 进入`example01_kumar`文件夹,运行Shell测试脚本`SX_Deconv_example_kumar.sh`

```bash
cd ./example01_kumar
bash ./SX_Deconv_example_kumar.sh
cd ../
```

结果放置于`example01_kumar/figs`文件夹中,图件如下：   

![Example Kumar](./example01_kumar/figs/corrected-XR.ST09.1998.197.121411-SRF.l.jpg)

kumar et al., 2006 Figure 4 原图见 [kumar et al., 2006](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JB003930)

#### 案例测试02： 基于脉冲反褶积的GC_SRF计算流程

- 进入`example02_gcsrf`文件夹,运行Shell测试脚本`GC_SRF_example_syn.sh`

该脚本中使用`GC_SRF`和`GC_SRF_post`程序提供了基于`GC_SRF`开展S波接收函数计算的参考流程.

```bash
cd ./example02_gcsrf
bash ./GC_SRF_example_syn.sh
```   

- 然后使用python执行绘图脚本

```bash
python ./GC_SRF_plot.py
```   

结果图件放置于`example02_gcsrf/figs`文件夹中,图件如下：   

![GC_SRF 001](./example02_gcsrf/figs/syn_SM_d75_001.jpg)   
![GC_SRF 002](./example02_gcsrf/figs/syn_SM_d75_002.jpg)

#### 必要工具和库

##### Linux

1. Taup
2. GMT 6
3. SAC
4. MPIC++
5. GCC

##### Python Lib

1. obspy 
2. pandas 
3. numpy 
4. matplotlib

---

**注意**：此双语 README 由 Grok 3 生成。
