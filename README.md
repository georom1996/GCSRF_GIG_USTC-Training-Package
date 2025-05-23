# <img src="./logo.png" width="100"/>  GCSRF_GIG_USTC-Training-Package   

[中文](#中文版本) | [English](#english-version)

---

## 中文版本

**2023-地震学算法与程序培训班-S波接收函数提取策略GC_SRF相关文件**

### 📁 文件树表

**1. SX_Deconv和GC_SRF程序及案例**
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

**2. 单台S波接收函数的RMSE排序和筛选脚本及案例**
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

**3. 基于PSDM的S波接收函数的偏移成像python执行脚本及案例**
- Section003/  
    - README.md
    - example001_pyPSDM_runner_CDMOD.py
    - example002_pyPSDM_runner_IASP.py
    - data_gcsrf/
    - figs/
    - mod/
    - psdm/

### 🚀 使用方法

每个Section中分别的README.md文件中提供了具体依赖需求和案例操作流程，同行可以根据实际情况选择使用。即便如此，仍然建议在Anaconda/Miniconda下安装大部分依赖程序，为方便各位使用，该部分将尽可能给出一个能够满足本程序包完整运行的最小的依赖安装方案。设备中若已有conda+sac则大概率能完整运行。

#### 1. 安装 conda 和 SAC

**conda**安装包链接见[conda website](https://docs.conda.io/en/latest/miniconda.html)

例如，如果选择Miniconda3 Linux 64-bit平台，则参考如下命令安装：
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

安装过程中全程使用默认选项即可，安装结束后需要初始化一下环境变量使conda生效，以bash脚本配置为例，参考如下命令初始化：
```bash
source ~/.bashrc
```

**SAC**的安装流程参考[SeismanBlog](https://seisman.github.io/SAC_Docs_zh/install/)。

#### 2. 配置 conda 和所需编译环境

**方案一：使用GCSRF.yml文件**
1. 下载最新版本中根目录中提供的GCSRF.yml文件，在程序包根目录中进行环境配置
2. 在conda中创建GCSRF环境并安装相应依赖：
   ```bash
   conda env create -f GCSRF.yml
   ```
3. 安装Taup：
   ```bash
   sudo apt-get install default-jdk
   sudo snap install taup
   ```

**方案二：手动安装依赖**
1. 可在自己常用的conda环境中安装如下依赖：
   ```bash
   conda config --add channels conda-forge
   conda install obspy=1.4.0
   conda install pandas
   conda install octave -c conda-forge
   conda install gmt -c conda-forge
   ```
2. 安装Taup：
   ```bash
   sudo apt-get install default-jdk
   sudo snap install taup
   ```

**编译环境检查：**
若Section001/src/中编译顺利，请跳过此部分。

检查编译器是否存在：
```bash
g++ --version
gcc --version
gfortran --version
mpic++ --version
```

若命令不存在则需根据需要安装：
- **Ubuntu或其他基于Debian的Linux发行版：**
  ```bash
  sudo apt update
  sudo apt install gcc gfortran mpich
  ```
- **Fedora, CentOS, RHEL等基于RPM的Linux发行版：**
  ```bash
  sudo dnf install gcc gfortran mpich
  ```
- **macOS：**
  ```bash
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
  brew install gcc gfortran mpich
  ```

### ❓ 常见问题 & 解决方案

<details>
<summary>点击展开查看常见问题</summary>

**问题1：** 编译时屏幕显示 `sac_lpz.c: 492:2 for loop initial declaration are only allowed in C99 mode In function 'read_data_':`

**解决方案：** 可在Makefile中为cpp文件的编译命令添加 -std=c99

**问题2：** `example01_kumar.sh` 运行时出现calculate文件夹不存在的提示

**解决方案：** 脚本中已更改并加入文件夹判断，错误已更正

**问题3：** GlibC和libstdc++版本问题可能导致在Redhat/CentOS平台编译失败

**解决方案：** 对于Ubuntu/Debian平台，可以通过更新libc版本解决：
```bash
sudo apt-get update
sudo apt-get install libc6
```

**更多问题请参考完整文档...**

</details>

### 📚 参考文献

**GC_SRF策略，请引用：**
- Zhang Z., & Deng Y. (2022). A Generalized Strategy From S‐Wave Receiver Functions Reveals Distinct Lateral Variations of Lithospheric Thickness in Southeastern Tibet. Geochemistry, Geophysics, Geosystems, 23(11). https://doi.org/10.1029/2022GC010619

**脉冲反褶积计算，请附加引用：**
- Treitel, S., & Robinson, E. (1966). Seismic wave propagation in layered media in terms of communication theory. Geophysics, 31(1), 17–32.
- Robinson, E., & Treitel, S. (1976). Net downgoing energy and resulting minimum-phase property of downgoing waves. Geophysics, 41(6), 1394–1396.
- Robinson, E., & Treitel, S. (2000). Geophysical signal analysis. Society of Exploration Geophysicists.

**PSDM成像，请附加引用：**
- Chen, L., Wen, L., & Zheng, T. (2005a). A wave equation migration method for receiver function imaging: 1. Theory. Journal of Geophysical Research, 110(B11).
- Chen, L., Wen, L., & Zheng, T. (2005b). A wave equation migration method for receiver function imaging: 2. Application to the Japan subduction zone. Journal of Geophysical Research, 110(B11).

### 🙏 致谢

现有GC_SRF流程的相关功能实现均有赖于同行学者的慷慨分享与交流讨论。特别感谢刘启民、徐强、Lupei Zhu、LP Vinnik、Xiaohui Yuan、陈凌研究员、吴庆举研究员、姚华建教授、Rongjiang Wang教授等学者的帮助，以及GIG Geophysics Lab - 724A团队和众多同行的支持。

### 🔮 展望

当前版本基于实际工作流程进行整理，暂未对脚本语言进行完全地统一，相关工作已在开展。目前各项功能函数均已在python环境中测试实现，更便于使用的新版本请关注后续PyGCSRF的上线。

---

## English Version

**2023 Seismological Algorithm and Programming Training - GC_SRF S-wave Receiver Function Extraction Strategy Related Files**

### 📁 File Structure

**1. SX_Deconv and GC_SRF Programs with Examples**
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

**2. Single Station S-wave Receiver Function RMSE Sorting and Filtering Scripts with Examples**
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

**3. PSDM-based S-wave Receiver Function Migration Imaging Python Scripts with Examples**
- Section003/  
    - README.md
    - example001_pyPSDM_runner_CDMOD.py
    - example002_pyPSDM_runner_IASP.py
    - data_gcsrf/
    - figs/
    - mod/
    - psdm/

### 🚀 Usage Instructions

Each Section contains detailed dependency requirements and case operation procedures in separate README.md files. Users can choose to use them according to their actual situation. Nevertheless, it is still recommended to install most dependent programs under Anaconda/Miniconda. For convenience, this section will provide a minimal dependency installation scheme that can satisfy the complete operation of this package. If the device already has conda+sac, it can likely run completely.

#### 1. Install conda and SAC

**conda** installation package link: [conda website](https://docs.conda.io/en/latest/miniconda.html)

For example, if choosing Miniconda3 Linux 64-bit platform, refer to the following commands for installation:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Use default options throughout the installation process. After installation, initialize environment variables to make conda effective. For bash script configuration, initialize with:
```bash
source ~/.bashrc
```

**SAC** installation process reference: [SeismanBlog](https://seisman.github.io/SAC_Docs_zh/install/).

#### 2. Configure conda and Required Compilation Environment

**Option 1: Using GCSRF.yml file**
1. Download the GCSRF.yml file provided in the root directory of the latest version
2. Create GCSRF environment in conda and install dependencies:
   ```bash
   conda env create -f GCSRF.yml
   ```
3. Install Taup:
   ```bash
   sudo apt-get install default-jdk
   sudo snap install taup
   ```

**Option 2: Manual dependency installation**
1. Install the following dependencies in your commonly used conda environment:
   ```bash
   conda config --add channels conda-forge
   conda install obspy=1.4.0
   conda install pandas
   conda install octave -c conda-forge
   conda install gmt -c conda-forge
   ```
2. Install Taup:
   ```bash
   sudo apt-get install default-jdk
   sudo snap install taup
   ```

**Compilation Environment Check:**
If compilation in Section001/src/ is successful, please skip this part.

Check if compilers exist:
```bash
g++ --version
gcc --version
gfortran --version
mpic++ --version
```

If commands don't exist, install as needed:
- **Ubuntu or other Debian-based Linux distributions:**
  ```bash
  sudo apt update
  sudo apt install gcc gfortran mpich
  ```
- **Fedora, CentOS, RHEL and other RPM-based Linux distributions:**
  ```bash
  sudo dnf install gcc gfortran mpich
  ```
- **macOS:**
  ```bash
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
  brew install gcc gfortran mpich
  ```

### ❓ Common Issues & Solutions

<details>
<summary>Click to expand common issues</summary>

**Issue 1:** Compilation shows `sac_lpz.c: 492:2 for loop initial declaration are only allowed in C99 mode In function 'read_data_':`

**Solution:** Add -std=c99 to the cpp file compilation command in Makefile

**Issue 2:** `example01_kumar.sh` shows calculate folder not found

**Solution:** Script has been modified with folder checking, error corrected

**Issue 3:** GlibC and libstdc++ version issues may cause compilation failure on Redhat/CentOS platforms

**Solution:** For Ubuntu/Debian platforms, resolve by updating libc version:
```bash
sudo apt-get update
sudo apt-get install libc6
```

**For more issues, please refer to the complete documentation...**

</details>

### 📚 References

**For GC_SRF strategy, please cite:**
- Zhang Z., & Deng Y. (2022). A Generalized Strategy From S‐Wave Receiver Functions Reveals Distinct Lateral Variations of Lithospheric Thickness in Southeastern Tibet. Geochemistry, Geophysics, Geosystems, 23(11). https://doi.org/10.1029/2022GC010619

**For spike deconvolution calculation, additionally cite:**
- Treitel, S., & Robinson, E. (1966). Seismic wave propagation in layered media in terms of communication theory. Geophysics, 31(1), 17–32.
- Robinson, E., & Treitel, S. (1976). Net downgoing energy and resulting minimum-phase property of downgoing waves. Geophysics, 41(6), 1394–1396.
- Robinson, E., & Treitel, S. (2000). Geophysical signal analysis. Society of Exploration Geophysicists.

**For PSDM imaging, additionally cite:**
- Chen, L., Wen, L., & Zheng, T. (2005a). A wave equation migration method for receiver function imaging: 1. Theory. Journal of Geophysical Research, 110(B11).
- Chen, L., Wen, L., & Zheng, T. (2005b). A wave equation migration method for receiver function imaging: 2. Application to the Japan subduction zone. Journal of Geophysical Research, 110(B11).

### 🙏 Acknowledgments

The implementation of GC_SRF workflow functions relies on the generous sharing and discussion exchanges of fellow scholars. Special thanks to Liu Qimin, Xu Qiang, Lupei Zhu, LP Vinnik, Xiaohui Yuan, Researcher Chen Ling, Researcher Wu Qingju, Professor Yao Huajian, Professor Rongjiang Wang, and others for their help, as well as the support from the GIG Geophysics Lab - 724A team and many colleagues.

### 🔮 Future Prospects

The current version is organized based on actual workflow procedures, so script languages have not been completely unified yet. Related work is underway. Currently, all functional functions have been tested and implemented in the Python environment. Please stay tuned for the upcoming PyGCSRF for a more user-friendly new version.
