## 操作流程
- 进入Section002/data/文件夹中，对单台数据进行RMSE排序操作，linux执行命令
```bash
cd ./data/
bash Do_RMSE.sh
cd ../
```
该脚本中使用了octave的静默运行，需提前安装

Ubuntu/WSL平台安装命令为```sudo apt install octave-signal```  
Anaconda/miniconda平台安装命令为```conda install -c conda-forge octave ```  
macOS Homebrew安装命令为```brew install octave```  
macOS App Bundles安装包见 https://github.com/octave-app/octave-app/releases


- 在Section002文件夹中，linux执行命令   
```bash
bash ./Plot_RMSE.sh CQ.ROC
```
结果图件放置于figs文件夹中，图件如下：   
![RMSE 001](./figs/CQ.ROC.Fig07.rmse.jpg)     

- 在Section002文件夹中，linux执行命令   
```bash
bash ./Plot_RMSE.sh SC.MGU
```
结果图件放置于figs文件夹中，图件如下：
![RMSE 002](./figs/SC.MGU.Fig07.rmse.jpg)   

## 必要工具和库
### Linux
1. GMT 6
2. SAC
3. Taup
4. octave