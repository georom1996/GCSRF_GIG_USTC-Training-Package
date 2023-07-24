# <img src="./logo.png" width="100"/>  GCSRF_GIG_USTC-Training-Package   

2023-地震学算法与程序培训班-S波接收函数提取策略GC_SRF相关文件   

### 文件树表 ：

***1. SX_Deconv和GC_SRF程序及案例***
- Section001/  
    - README.md
    - src/
    - bin/
    - example01_kumar/
        - SX_Deconv_example_kumar.sh
        - calculate_test/
        - data/
        - figs/
    - example02_gcsrf/
        - GC_SRF_example_syn.sh
        - GC_SRF_plot.py
        - calculate_test/
        - data/
        - figs/
        - gcsrf_plot_subroutines.py   

***2.单台S波接收函数的RMSE排序和筛选脚本及案例***
- Section002/  
    - Plot_RMSE.sh
    - README.md
    - data/
        - CQ.ROC/
        - SC.MGU/
        - Do_RMSE.sh
        - fdd.lst
        - readsac.m
        - Silent_FiltRMSE.m

    - figs/
    - tmp/   

***3. 基于PSDM的S波接收函数的偏移成像python执行脚本及案例***
- Section003/  
    - README.md
    - example001_pyPSDM_runner_CDMOD.py
    - example002_pyPSDM_runner_IASP.py
    - data_gcsrf/
    - figs/
    - mod/
    - psdm/


### 使用方法 ：

- 每个Section中分别的README.md文件中提供了具体依赖需求和案例操作流程，同行可以根据实际情况选择使用。
- 即便如此，仍然建议在Anaconda/Miniconda下安装大部分依赖程序，为方便各位使用，该部分将尽可能给出一个能够满足本程序包完整运行的最小的依赖安装方案。
- 设备中若已有conda+sac则大概率能完整运行。

**conda 方案一：**
1. 下载最新版本中根目录中提供的GCSRF.yml文件
2. 在conda中创建GCSRF环境并安装相应依赖：
    ```
    conda env create -f GCSRF.yml
    ```
3. 安装Taup:
    ```
    sudo apt-get install default-jdk
    sudo snap install taup
    ```

**conda 方案二：**
1. 在自己常用的conda环境中安装如下依赖：
    ```
    conda config --add channels conda-forge
    conda install obspy=1.4.0
    conda install pandas
    conda install octave -c conda-forge
    conda install gmt -c conda-forge
    ```
3. 安装Taup:
    ```
    sudo apt-get install default-jdk
    sudo snap install taup
    ```
**编译工具方案：**
- 除这些之外，还需检查gcc, gfortran, g++, mpich等编译器是否存在:
    ```
    g++ --version
    gcc --version
    gfortran --version
    mpic++ --version
    ```

    若命令不存在则需根据需要安装:
    **> Ubuntu或其他基于Debian的Linux发行版:**
    打开终端并输入以下命令：
    ```
    sudo apt update
    sudo apt install g++
    ```
    **> Fedora, CentOS, RHEL等基于RPM的Linux发行版：**

    打开终端并输入以下命令：
    ```
    sudo dnf install gcc-c++
    ```
    **> macOS：**
    需安装Homebrew后使用brew命令安装：
    ```
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
    brew install gcc
    ```


### 常见问题 & 解决方案 2023-07-23

**I1:**

- ```编译时屏幕显示 sac_lpz.c: 492:2  for loop initial declaration are only allowed in C99 mode In function 'read_data_':``` 

**S1:**

- 可在Makefile 中为cpp文件的编译命令添加 -std=c99 **【Ubuntu测试无需，其他平台根据需要使用，最新Makefile中提供了这一选项】**

**I2:**

- ```example01_kumar.sh 运行时出现calculate文件夹不存在的提示```

**S2:**

- 脚本中已更改并加入文件夹判断，错误已更正。 **【已修改】**

**I3:**

- ```GlibC和libstdc++版本问题可能导致在Redhat/CentOS平台编译失败```

**S3:**

- 对于Ubuntu/Debian平台，可以通过更新libc版本解决：

    ```
    sudo apt-get update
    sudo apt-get install libc6
    ```

**I4:**

- ```编译提示libgfortran.so.5找不到```

**S4:**
- 对于Ubuntu/Debian平台，可以通过安装libgfortran包解决：

    ```
    sudo apt-get update
    sudo apt-get install libgfortran5
    ```


**I5:**
- 该问题详见 [Closed Issue1](https://github.com/georom1996/GCSRF_GIG_USTC-Training-Package/issues/1)
    ```
    PSDM源代码中可能存在数组越界的问题：
    There might be a possible out-of-range of the array **inn1** in pierc_new_n.f: 297 line.
    ```

**S5:**

- 为避免出现数组越界问题，已对相应变量进行了初始化。**【已修改】**

**I6:**

- ```执行Section001中第一个例子时出现pssac WARNING，显示tmarker9没有定义```

**S6:**
- 原因在于shell中未能正确调用taup_setsac。若Taup是使用```sudo snap install taup```安装，则需确保java安装正常，若出现java库不存在的提示，可尝试通过```sudo apt-get install default-jdk```解决。

**I7:**

- ```在部分平台编译时需要指定gcc版本后才能编译，，f77可以替换成gfortran```

**S7:**

- 测试显示g++-11和g++-10均能编译通过 **【f77已替换为gfortran】**

**I8:**

- ```“omp.h” 找不到```

**S8:**

- 程序的早期版本使用过omp，当前版本未使用该库。**【已注释该行】**



### 参考文献 ：
GC_SRF策略，请引用：   
- Zhang Z., & Deng Y. (2022). A Generalized Strategy From S‐Wave Receiver Functions Reveals Distinct Lateral Variations of Lithospheric Thickness in Southeastern Tibet. Geochemistry, Geophysics, Geosystems, 23(11). https://doi.org/10.1029/2022GC010619


脉冲反褶积计算，请附加引用：   
- Treitel, S., & Robinson, E. (1966). Seismic wave propagation in layered media in terms of communication theory. Geophysics, 31(1), 17– 32. https://doi.org/10.1190/1.1439729   
- Robinson, E., & Treitel, S. (1976). Net downgoing energy and resulting minimum-phase property of downgoing waves. Geophysics, 41(6), 1394– 1396. https://doi.org/10.1190/1.1440689   
- Robinson, E., & Treitel, S. (2000). Geophysical signal analysis. Society of Exploration Geophysicists. https://doi.org/10.1190/1.9781560802327   

PSDM成像，请附加引用：   

- Chen, L., Wen, L., & Zheng, T. (2005a). A wave equation migration method for receiver function imaging: 1. Theory. Journal of Geophysical Research, 110(B11). https://doi.org/10.1029/2005JB003665   
- Chen, L., Wen, L., & Zheng, T. (2005b). A wave equation migration method for receiver function imaging: 2. Application to the Japan subduction zone. Journal of Geophysical Research, 110(B11). https://doi.org/10.1029/2005JB003666   

### 致谢 ：   
现有GC_SRF流程的相关功能实现均有赖于同行学者的慷慨分享与交流讨论。   

其中：   

- SX_Deconv和GC_SRF源码编写的早期工作得到了刘启民，徐强，Lupei Zhu等学者的帮助；   
- 接收函数基础以及QSEIS模型参数设置等方面得到了 LP Vinnik，Xiaohui Yuan 等学者的帮助；   
- 该工作的雏形也曾在2019年地震学算法与程序培训班(大理)上得到过陈凌研究员，吴庆举研究员，姚华建教授的当面指导；   
- QSEIS波形模拟程序来自 Rongjiang Wang 教授课题组;   
- PSDM相关程序来自于陈凌研究员课题组;   
- sacsnr、sacStack、sacio等均来自 Lupei Zhu 教授课题组； 
- 程序整理得到了 GIG Geophysics Lab - 724A 成员 ( 胡仲发、朱晟、黄润青、罗恒 ) 的帮助；
- 程序问题和解决方案部分得到了肖卓、韩如冰、苏文君柳、曾彦迪、闫晓东、寇华东等人的帮助。

在此一并表示感谢。


### 展望 ：
- 当前版本基于实际工作流程进行整理，故暂未对脚本语言进行完全地统一，相关工作已在开展。   
- 目前各项功能函数均已在python环境中测试实现，更便于使用的新版本请关注后续PyGCSRF的上线。