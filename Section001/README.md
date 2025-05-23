# Pulse Deconvolution and GC_SRF Strategy Implementation

<div style="text-align: center; margin-bottom: 20px;">
  <button onclick="showLanguage('en')" style="padding: 10px 20px; margin-right: 10px;">English</button>
  <button onclick="showLanguage('zh')" style="padding: 10px 20px;">中文</button>
</div>

<div id="en-content" style="display: block;">
## Compile the Base Program
- Navigate to the `src` folder and compile the source code:
```bash
cd src
make clean
make all
cd ../
```

### Executable Programs Description
After compiling, the following executable programs will be available in the `bin` folder:

- **1. SX_Deconv**
  This program is an independent computation tool for pulse deconvolution. Usage:
  ```bash
  SX_Deconv -Z Z.sac -N N.sac -E E.sac –inc 16 -t1 100 -t2 120 -D
  ```
  Program options:
  ```
  -t1 Pulse start time (relative to the SAC file start time)
  -t2 Pulse end time (relative to the SAC file start time)
  -inc Rotation angle
  -D Output inverse filter and rotated waveform SAC files
  ```
  Input files:
  ```
  Z.sac N.sac E.sac
  Note: Ensure the header contains baz
  ```
  Output files:
  ```
  Default output:
    rf.sl.sac : Non-normalized L-component deconvolution waveform (SRF)
    rf.sq.sac : Non-normalized Q-component deconvolution waveform
    rf.st.sac : Non-normalized T-component deconvolution waveform
    rf.nsl.sac : Amplitude-normalized L-component deconvolution waveform (SRF)
    rf.nsq.sac : Amplitude-normalized Q-component deconvolution waveform
    rf.nst.sac : Amplitude-normalized T-component deconvolution waveform
  If -D option is used, additional outputs:
    l-com.sac : L-component original waveform
    q-com.sac : Q-component original waveform
    t-com.sac : T-component original waveform
    fil.sac : Inverse wavelet filter
  ```

- **2. GC_SRF**
  This program iteratively computes deconvolution for different wavelet window lengths and rotation angles. Usage:
  ```bash
  GC_SRF -Z Z.sac -N N.sac -E E.sac -t1 t1
  ```
  Program options:
  ```
  -t1 Pulse start time (relative to the SAC file start time)
  Note: The program defaults to looping over win_len=0:5:100; inci=0:4:60, where t1 is the direct S arrival time relative to the SAC file start, and win_len starts at t1-5s
  ```
  Input files:
  ```
  Z.sac N.sac E.sac
  Note: Ensure the header contains baz
  ```
  Output files:
  ```
  Default output:
    ("rf-%02.f-%03.f.sl.sac", incnum, win_len) : Non-normalized L-component deconvolution waveform (SRF)
    ("rf-%02.f-%03.f.sq.sac", incnum, win_len) : Non-normalized Q-component deconvolution waveform
    ("rf-%02.f-%03.f.st.sac", incnum, win_len) : Non-normalized T-component deconvolution waveform
    ("rf-%02.f-%03.f.nsl.sac", incnum, win_len) : Amplitude-normalized L-component deconvolution waveform (SRF)
    ("rf-%02.f-%03.f.nsq.sac", incnum, win_len) : Amplitude-normalized Q-component deconvolution waveform
    ("rf-%02.f-%03.f.nst.sac", incnum, win_len) : Amplitude-normalized T-component deconvolution waveform
  ```

- **3. GC_SRF_par**
  This is the parallel version of GC_SRF, with the same usage, designed to improve the efficiency of iterative computations:
  ```bash
  GC_SRF_par -Z Z.sac -N N.sac -E E.sac -t1 t1
  ```
  Program options:
  ```
  -t1 Pulse start time (relative to the SAC file start time)
  -D Output inverse filter and rotated waveform SAC files
  ```

- **4. GC_SRF_post**
  This post-processing program for GC_SRF/GC_SRF_par computes the amplitude and energy grids at time 0 and cross-correlation coefficient grids for all deconvolution results:
  ```bash
  GC_SRF_post
  ```
  Program options:
  ```
  No options; run directly in the GC_SRF computation directory
  ```
  Input files:
  ```
  SAC format data files ending with nsl.sac.100.cutrv.rm.rv, processed from GC_SRF/GC_SRF_par output
  Note: To obtain these input SAC files:
  1. Truncate GC_SRF/GC_SRF_par output waveforms (nsl.sac files) from -100s to 100s (.100.cut)
  2. Reverse the time axis (rv)
  3. Apply moveout correction (.rm)
  4. Reverse the time axis again (rv)
  Refer to the script in Section001/example02_gcsrf/GC_SRF_example_syn.sh
  ```
  Output files:
  ```
  all.Amp_0.xyz : Amplitude spectrum of all SRFs at time 0
  all.Energy_0.xyz : Energy spectrum of all SRFs at time 0
  Min.lst : List of waveforms with the lowest amplitude at time 0 for different window lengths, format: rotation angle, window length, filename
  all.CCC.xyz : Cross-correlation coefficient spectrum of all SRFs with the mean waveform
  SelectedRed.CC.lst : Cross-correlation coefficients for waveforms with the lowest amplitude at time 0, format: rotation angle, window length, filename
  OptimalSRF.lst : Optimal SRF, format: filename, rotation angle, window length, cross-correlation coefficient
  ```

- **5. spmove.joe**
  ```
  Modified from Professor Xiaohui Yuan's spmout code
  ```

- **6. sacsnr and sacStack**
  From [SACTOOLS](https://github.com/msthorne/SACTOOLS/tree/master).

### Case Study 01: Reproducing Kumar et al., 2006 Figure 4b
- Navigate to the `example01_kumar` folder and run the shell test script `SX_Deconv_example_kumar.sh`:
```bash
cd ./example01_kumar
bash ./SX_Deconv_example_kumar.sh
cd ../
```
Results are stored in the `example01_kumar/figs` folder, with the following figure:
![Example Kumar](./example01_kumar/figs/corrected-XR.ST09.1998.197.121411-SRF.l.jpg)

Original figure from Kumar et al., 2006: [Kumar et al., 2006](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JB003930)

### Case Study 02: GC_SRF Computation Workflow Based on Pulse Deconvolution
- Navigate to the `example02_gcsrf` folder and run the shell test script `GC_SRF_example_syn.sh`. This script uses GC_SRF and GC_SRF_post to provide a reference workflow for S-wave receiver function computation:
```bash
cd ./example02_gcsrf
bash ./GC_SRF_example_syn.sh
```
- Then, execute the plotting script using Python:
```bash
python ./GC_SRF_plot.py
```
Result figures are stored in the `example02_gcsrf/figs` folder, with the following figures:
![GC_SRF 001](./example02_gcsrf/figs/syn_SM_d75_001.jpg)
![GC_SRF 002](./example02_gcsrf/figs/syn_SM_d75_002.jpg)

## Required Tools and Libraries
### Linux
1. Taup
2. GMT 6
3. SAC
4. MPIC++
5. GCC

### Python Libraries
1. obspy
2. pandas
3. numpy
4. matplotlib

</div>

<div id="zh-content" style="display: none;">
## Section001 脉冲反褶积与GC_SRF策略的实现
## 基础程序编译
- 进入src文件夹,编译相关源码
```bash
cd src
make clean
make all
cd ../
```

### 可执行程序说明
在流程操作编译完成之后,在bin文件夹中将出现如下可执行程序：
- **1. SX_Deconv**
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

- **2. GC_SRF**
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

- **3. GC_SRF_par**
  该程序为GC_SRF的并行版本,使用方法与GC_SRF相同,用于加快循环计算的效率：
  ```bash
  GC_SRF_par -Z Z.sac -N N.sac -E E.sac -t1 t1
  ```
  程序选项：
  ```
  -t1 脉冲起始时刻 (相对SAC文件开始时刻的时间)
  -D 输出逆滤波器及波形旋转后的SAC文件
  ```

- **4. GC_SRF_post**
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

- **5. spmove.joe**
  ```
  修改自Xiaohui Yuan教授spmout代码
  ```

- **6. sacsnr 和 sacStack**
  来自[SACTOOLS](https://github.com/msthorne/SACTOOLS/tree/master).

### 案例测试01：复现Kumar et al., 2006 Figure 4b
- 进入example01_kumar文件夹,运行Shell测试脚本SX_Deconv_example_kumar.sh
```bash
cd ./example01_kumar
bash ./SX_Deconv_example_kumar.sh
cd ../
```
结果放置于example01_kumar/figs文件夹中,图件如下：
![Example Kumar](./example01_kumar/figs/corrected-XR.ST09.1998.197.121411-SRF.l.jpg)

Kumar et al., 2006 Figure 4 原图见 [Kumar et al., 2006](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JB003930)

### 案例测试02：基于脉冲反褶积的GC_SRF计算流程
- 进入example02_gcsrf文件夹,运行Shell测试脚本GC_SRF_example_syn.sh
该脚本中使用GC_SRF和GC_SRF_post程序提供了基于GC_SRF开展S波接收函数计算的参考流程.
```bash
cd ./example02_gcsrf
bash ./GC_SRF_example_syn.sh
```
- 然后使用python执行绘图脚本
```bash
python ./GC_SRF_plot.py
```
结果图件放置于example02_gcsrf/figs文件夹中,图件如下：
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

</div>

## Statement
This bilingual README was generated by Grok 3.

<script>
function showLanguage(lang) {
  document.getElementById('en-content').style.display = lang === 'en' ? 'block' : 'none';
  document.getElementById('zh-content').style.display = lang === 'zh' ? 'block' : 'none';
}
</script>
