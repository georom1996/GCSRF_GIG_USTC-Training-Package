# README

本仓库提供中英文两种语言的文档说明。

- [English](#english-version)
- [Chinese](#chinese-version)

## English Version <a name="english-version"></a>

### Section003 PSDM Offset Imaging of Multi-Station S-Wave Receiver Functions

#### Script Explanation

**Function:**

- A Python Runner built based on the PSDM program package from Chen Ling's research group. In this Runner script, the parameter settings are consistent with those of the PSDM program package. The parameter setting section is identified by the fields "# [Parameters Part Begin]" and "# [Parameters Part End]".

#### Offset Imaging Operation Flow

- In the `Section003` folder, run the Linux command:
  ```bash
  python example001_pyPSDM_runner_CDMOD.py
  ```
  The resulting figures are placed in the `figs` folder, as shown:
  ![PSDM CDMOD](./figs/example001_CDmod_PSDM_Map.jpg)

- In the `Section003` folder, run the Linux command:
  ```bash
  python example002_pyPSDM_runner_IASP.py
  ```
  The resulting figures are placed in the `figs` folder, as shown:
  ![PSDM IASP](./figs/example002_IASP_PSDM_Map.jpg)

#### Necessary Tools and Libraries

##### Linux
1. SAC

##### Python
1. obspy
2. numpy
3. pandas
4. scipy
5. matplotlib

## Chinese Version <a name="chinese-version"></a>

### Section003 多台S波接收函数的PSDM偏移成像

#### 脚本说明

**功能:**

- 以陈凌老师课题组PSDM程序包为基础构建的python Runner,在该Runner脚本中相应参数设置与PSDM程序包一致.参数设置部分由"# [Parameters Part Begin]"和"# [Parameters Part End]"字段标识.

#### 偏移成像操作流程

- 在Section003文件夹中，linux执行命令   
```bash
python example001_pyPSDM_runner_CDMOD.py
```
结果图件放置于figs文件夹中，图件如下:   
  ![PSDM IASP](./figs/example002_IASP_PSDM_Map.jpg)
