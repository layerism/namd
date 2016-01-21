# namd

## 新增加了ITS, ACG 采样模块

**积分回火增强采样算法(ITS)**: [An integrate-over-temperature approach for enhanced sampling] (http://scitation.aip.org/content/aip/journal/jcp/128/6/10.1063/1.2825614)

修改的部分用均用

//lyk

... ...

//lyk 

包围

## ITS模块配置文件参数解释

```sh
its                   on		    # on/off, 激活/关闭ITS算法
itsFile               its.pdb   # pdb文件, 用于指定原子分类
itsCol                B         # 指定某列用于原子分类, 比如B因子列
groupIdx1             0         # 类1的指标, 原子分类类别1的类标
groupIdx2             2         # 类2的指标, 原子分类类别2的类
itsOutFreq            1000      # ITS输出间隔, 分子模拟步数
itsFirstStep          0         # 何时开始ITS
itsLastStep           100000000 # 合适结束ITS
itsWhamTimes          200       # WHAM的更新次数
itsUpdateFreq         100000    # 每次WHAM的运算步长
itsSamplingFreq       25        # 采样间隔

itsUpdate             on        # 是否自动更新ITS参数
itsUpdateLog          ITSUPDATELOG.DAT  # 更新参数日记
#itsInFile            ITSINPUT          # 如果指明, ITS参数初始值将从这里读入, 不然便采用指定范围内的均匀分布
itsReplicaNum         21	              # ITS副本个数
itsLbdMin             -0.5              # ITS温度因子重标的最小值
itsLbdMax             0.0               # ITS温度因子重标的最大值
itsIntType            vdw               # ITS算法只作用于vdw相互作用, 还有elect, elect-vdw选项, 表示短程电, 短程电+vdw
onlyInt               off               # on表示ITS只用到类1和类2的相互作用上, 自身与自身的作用部分不执行ITS
itsTemp               300               # ITS的工作温度, 通常取同系统模拟温度一样		            
						
DENSITYALAPHA         2.0
WCAP                  5.0
WHAMSKIP              1
WHAMLOGZERROR         2.0
```



