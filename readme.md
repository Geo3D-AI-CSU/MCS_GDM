
[Uncertainty Evaluation of Soil Heavy Metal(loid) Pollution and Health Risk: A Geographic Detector with Monte Carlo Simulation](https://www.mdpi.com/2305-6304/11/12/1006 "论文")

# 环境
python==3.7  
pymc== 2.3.8  
numpy==1.21.6  
gdal==2.2.4  

# 软件
arcmap

# 数据
Raw shp 包含原始土壤采样点和湖南省边界数据  
influence factors 包含16个影响因子栅格数据  
Monte Carlo output 包含使用MC模拟的采样点数据shp、属性数据Trace_npy、影响因子栅格数据Raster5000、健康风险指数涉及的变量HQ_CR  
Nemerow1000 包含由模拟的采样点和属性计算的Nemerow指数，插值、重分类、计算信息熵的栅格数据  
geodector 包含GDM的输入和输出数据。  

# 程序入口
python main.py  
GDM代码：batch_derectorNemerow.R



