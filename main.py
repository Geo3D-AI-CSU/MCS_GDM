#!/usr/bin/env python
# _*_ coding:utf-8 _*_
__author__ = 'syc'


# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os
from osgeo import ogr, gdal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm
from osgeo import ogr
import GetRasterEnt as GRE
import PointMethod
import threading
import attributeSampling

global PATH,RAW_DATA,FACTORS_PATH,MC_PATH,NEMEROW_PATH,PCDW_2000,HUNAN,MC_SHP,MC_ATTRIBUTE,HQ_CR
PATH=r'E:\syc\MCS_GDM_Master'
RAW_DATA=PATH+r"\Raw shp"
FACTORS_PATH=PATH+r"\influence factors"
MC_PATH=PATH+r"\Monte Carlo output"
#Nemerow data
NEMEROW_PATH=PATH+r"\Nemerow1000"
#采样点数据
PCDW_2000=RAW_DATA+r'\12PCDW_2000.shp'
HUNAN=RAW_DATA+r"\hunan_prj.shp"
MC_SHP=MC_PATH+r"\shp"
#健康风险指数数据
MC_ATTRIBUTE=MC_PATH+r"\Trace_npy"
HQ_CR=MC_PATH+r"\HQ_CR"


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print("Hi, {0}".format(name))  # Press Ctrl+F8 to toggle the breakpoint.



#HQ、CR出图
def show_cdf(npy_value):
    #HQ_CR_npy
    HQ_CR= np.load(npy_value)
    FieldData=[]#5*4个属性字段
    for j in range(8,len(HQ_CR[0]),12):
        k=j
        for k in range(j,j+4):
            FeaData=[]#48811各要素
            for i in range(len(HQ_CR)):
                FeaData.append(HQ_CR[i][k]*10**3)
                #print(FeaData)
            FieldData.append(FeaData)
    FieldData=np.array(FieldData)

    #FieldData：['Cr','Cd','As','Pb','Hg']>['HI_A','HI_K','TRC_A','TRC_K']共20个字段
    #plot.subplots_adjust(left=0,top=1)
    for i0 in range(10,len(FieldData),4): 
        plot=plt.figure()
        ax1 = plot.add_subplot(1,1,1)   
        #计算各元素HI、CR时取消注释
        #TransField_A=FieldData[i0].reshape(-1)
        #计算THI时取消注释
        #TransField_A=FieldData[0].reshape(-1)+FieldData[4].reshape(-1)+FieldData[8].reshape(-1)+FieldData[12].reshape(-1)+FieldData[16].reshape(-1)
        #计算TCR时取消注释
        TransField_A=FieldData[2].reshape(-1)+FieldData[6].reshape(-1)+FieldData[10].reshape(-1)#TCR
        denominator_A=len(TransField_A)#分母数量
        Data_A=pd.Series(TransField_A)#将数据转换为Series利用分组频数计算
        Fre_A=Data_A.value_counts()
        Fre_sort_A=Fre_A.sort_index(axis=0,ascending=True)
        #print(Fre_sort_A)
        Fre_df_A=Fre_sort_A.reset_index()#将Series数据转换为DataFrame
        Fre_df_A[0]=Fre_df_A[0]/denominator_A#转换成概率
        Fre_df_A.columns=['Rds','Fre']
        #Cr HI(index 0):-2,Cr CR(index 2):1,Cd HI(4):-1, Cd CR(index 6):3,As HI(8):-3,As CR(10):2,Pb HI(12):-2,Pb 无TCR(14),Hg HI(16):-1,THI:-3,TCR:1
        Fre_df_A['Rds']=Fre_df_A['Rds']*10**1
        #Cd HQ(4):25,Cd TCR(index 6):10,Cr HQ(0):25,Cr CR(index 2):5,As HI(8):10,As TCR(10):40,Pb HQ:10,Hg HI(16):10,HI:10,TCR:10
        x_max=10
        #Cd HQ(4):1,Cd TCR(index 6):1,Cr HQ:1,Cr CR(index 2):0.4,As HI(8):0.2,As TCR(10):4,Pb HQ:0.3,Hg HI(16):0.4,HI:0.25,TCR:0.7
        a_Alocation=0.7
        #Cd HQ(4):10,Cd TCR(index 6):2,Cr HQ:8,Cr CR(index 2):0.8,As HI(8):3,As TCR(10):10,Pb (12)HQ:2.5,Hg HI(16):3,HI:4,TCR:3
        a_Klocation=3
        #Cd HQ(4):15,Cd TCR(index 6):6,Cr HQ:15,Cr CR(index 2):3,As HI(8):7,As TCR(10):25,Pb HQ:7,Hg HI(16):7,HI:7,TCR:6
        p_location=6
        a_A=np.mean(Fre_df_A['Rds'])

        plt.text(a_Alocation,0.6,round(a_A,2) , fontdict = {
                    'family': 'Calibri', # 标注文本字体
                    'fontsize': 12, # 文本大小
                    'color': 'blue',  # 文本颜色
                })
        Fre_df_A['cumsum']=np.cumsum(Fre_df_A['Fre'])  
        ax1.axvline(x=a_A,color="blue",linestyle ="--")
        ax1.plot(Fre_df_A['Rds'],Fre_df_A['cumsum'],label='Adult')
        #计算各元素HI、CR时取消注释        
        # TransField_K=FieldData[i0+1].reshape(-1)
        #计算THI时取消注释
        #TransField_K=FieldData[0+1].reshape(-1)+FieldData[4+1].reshape(-1)+FieldData[8+1].reshape(-1)+FieldData[12+1].reshape(-1)+FieldData[16+1].reshape(-1)
        #计算TCR时取消注释
        TransField_K=FieldData[2+1].reshape(-1)+FieldData[6+1].reshape(-1)+FieldData[10+1].reshape(-1)
        denominator_K=len(TransField_K)
        Data_K=pd.Series(TransField_K)
        Fre_K=Data_K.value_counts()
        Fre_sort_K=Fre_K.sort_index(axis=0,ascending=True)
        #print(Fre_sort_K)
        Fre_df_K=Fre_sort_K.reset_index()
        Fre_df_K[0]=Fre_df_K[0]/denominator_K
        #Cr HI(index 0):-2,Cr CR(index 2):1,Cd HI(4):-1,Cd CR(index 6):3,As HI(8):-3,As TCR(10):2,,Pb HI(12):-2,Hg HI(16):-1,THI:-3,TCR:1
        Fre_df_K['Rds']=Fre_df_K['Rds']*10**1
        a_K=np.mean(Fre_df_K['Rds'])
        print("mean_A:",a_A)
        print("mean_K:",a_K)
  
        plt.text(a_Klocation, 0.6,round(a_K,2) , fontdict = {
            'family': 'Calibri', # 标注文本字体
            'fontsize': 12, # 文本大小
            'color': 'orange',  # 文本颜色
        })
        Fre_df_K['cumsum']=np.cumsum(Fre_df_K['Fre'])
        ax1.axvline(x=a_K,color="orange",linestyle ="--")
        ax1.plot(Fre_df_K['Rds'],Fre_df_K['cumsum'],label='Kids')
        

        ax1.set_title("CDF")
        ax1.set_ylabel("P")
        ax1.set_xlim(0,x_max)
        
        ax1.legend(loc='upper right')            #显示图例,plt.legend()
        #每4个字段绘制一次HI
        if (i0//2%2)==0:
            print(Fre_df_A.loc[Fre_df_A['Rds'] < 1.0*10**0,'cumsum'])
            p=max(Fre_df_A.loc[Fre_df_A['Rds'] < 1.0*10**0,'cumsum'])#Cr HI(index 0):1,Cd HQ(4):2,Cd TCR(index 6):6,As HI(8):0,As TCR(10):5,Pb HI(12):1,Hg HI(16):2,HI:0
            p0=max(Fre_df_K.loc[Fre_df_K['Rds'] < 1.0*10**0,'cumsum'])
            plt.text(p_location, 0.2,'Adults: ({}%>1)'.format(round(1-p,2)*100) , fontdict = {
                    'family': 'Calibri', 
                    'fontsize': 12, 
                    'color': 'blue', 
                })
            plt.text(p_location, 0.1,'Kids: ({}%>1)'.format(round(1-p0,2)*100)  , fontdict = {
            'fontsize': 12,
            'color': 'orange', 
                })
            #Cr HI(index 0):-1,Cd HI(4):E-2,As HI(8):无  ,Pb HI(12):E-1,Hg HI(16):E-2,HI:0
            ax1.set_xlabel("HI")
            plt.show()
            plot.savefig(r'E:\syc\soil\HI_new\HI.png')
        else:
            #判断CR是否小于10**-4标准（10的-4需对应缩放）
            x0=Fre_df_A.loc[Fre_df_A['Rds'] < 1.0*10**0]
            #Cr CR(index 2):-2,Cd CR(index 6):0,As HI(8):-6,,As TCR(10):-1,Pb HI(12):-5,TCR:-2#10**-6标准
            x1=Fre_df_A.loc[Fre_df_A['Rds'] < 1.0*10**-2]
            if len(x1) != 0:
                p0=max(Fre_df_A.loc[Fre_df_A['Rds'] < 1.0*10**0,'cumsum'])
                print('Adults: ({}%>10E-4)'.format(round(1-p0,2)*100))
                p1=max(Fre_df_A.loc[Fre_df_A['Rds'] < 1.0*10**-2,'cumsum'])
                plt.text(p_location, 0.2,'Adults: ({}%>10E-6)'.format(round(1-p1,2)*100) , fontdict = {
                    'family': 'Calibri', 
                    'fontsize': 12,
                    'color': 'blue', 
                })
            #全部大于10-6的标准，则关注大于10-4的概率
            else:
                if len(x0)!=0: 
                    p0=max(Fre_df_A.loc[Fre_df_A['Rds'] < 1.0*10**0,'cumsum']) 
                    print('Adults: ({}%>10E-4)'.format(round(1-p0,2)*100))  
                plt.text(p_location, 0.2,'Adults: (100.00%>10E-6)', fontdict = {
                    'family': 'Calibri',
                    'fontsize': 12, 
                    'color': 'blue', 
                })
            
            x0=Fre_df_K.loc[Fre_df_K['Rds'] < 1.0*10**0]
            #Cr CR(index 2):-2,Cd CR(index 6):0,,As HI(8):-6,,As TCR(10):-1,,Pb HI(12):-5,TCR:-2
            x1=Fre_df_K.loc[Fre_df_K['Rds'] < 1.0*10**-2]
            if len(x1) != 0:          
                p0=max(Fre_df_K.loc[Fre_df_K['Rds'] < 1.0*10**0,'cumsum'])      
                p1=max(Fre_df_K.loc[Fre_df_K['Rds'] < 1.0*10**-2,'cumsum'])
                print('Kids: ({}%>10E-4)'.format(round(1-p0,2)*100))     
                plt.text(p_location, 0.1,'Kids: ({}%>10E-6)'.format(round(1-p1,2)*100)  , fontdict = {
                'family': 'Calibri', 
                'fontsize': 12, 
                'color': 'orange',
                    })
            else: 
                if len(x0)!=0: 
                    p0=max(Fre_df_K.loc[Fre_df_K['Rds'] < 1.0*10**0,'cumsum']) 
                    print('Kids: ({}%>10E-4)'.format(round(1-p0,2)*100))     
                plt.text(p_location, 0.1,'Kids: (100.00%>10E-6)' , fontdict = {
                'family': 'Calibri', 
                'fontsize': 12,
                'color': 'orange', 
                    })
            #Cr CR(index 2):E-4,Cd CR:E-6,As CR(10):E-5,TCR:E-4  每出一张图需对应更改为下一标注
            ax1.set_xlabel("TCR(E-4)")
            plt.show()
            plot.savefig(r'E:\syc\soil\HI_new\TCR.png')
        




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    
    FIELD_NAME_ARR =[u"土壤Cr",u"土壤Pb",u"土壤Cd",u"土壤As",u"土壤Hg",u"土壤pH"]
    #Raster data
    RASTER_NAME=['NemerowRec']

    # 生成采样点模拟数据
    PointMethod.Point_Method([PCDW_2000],["12PCDW_T"],MC_SHP)
    #生成属性模拟数据
    SIGMA_ARR = [0.2,0.2,0.2,0.2,0.2,0.02]
    attributeSampling.main(PCDW_2000,MC_ATTRIBUTE, FIELD_NAME_ARR, SIGMA_ARR)

    #计算内梅罗指数后
    attributeSampling.calNemerow(MC_SHP, '12PCDW_T','GB15816')

    #IDW插值生成内梅洛指数栅格图
    attributeSampling.IDWInter(MC_SHP+r"\12PCDW_T", NEMEROW_PATH, ["Nemerow"], ['Nemerow'])
    
    # 在插值结果进行重分类后计算信息熵
    GRE.EntMethod(NEMEROW_PATH,RASTER_NAME)
    #按掩膜提取Ent结果
    attributeSampling.ExtractMask(NEMEROW_PATH, HUNAN,NEMEROW_PATH, ["NemerowRec"])

    #模拟健康风险指数变量
    #BW
    BW_Anpy=HQ_CR+r"\BW_A.npy"
    BW_Knpy=HQ_CR+r"\BW_K.npy"
    if os.path.exists(BW_Anpy)==False:
        BW_ATrace=attributeSampling.LognormalSample(60.1,12.1)
        print(BW_ATrace)
        BW_KTrace=attributeSampling.LognormalSample(19.3,3.2)
        print(BW_KTrace)
        np.save(BW_Anpy, BW_ATrace)
        np.save(BW_Knpy, BW_KTrace)
    BW_ATrace = np.load(BW_Anpy)
    BW_KTrace = np.load(BW_Knpy)
    print('BW_A:')
    print("负数：",BW_ATrace[BW_ATrace<0])
    print("负数：",BW_KTrace[BW_KTrace<0])
    print('第50分位数:{}'.format(np.percentile(BW_ATrace,50)))
    print('第95分位数:{}'.format(np.percentile(BW_ATrace,95)))
    print('BW_K:')
    print('第50分位数:{}'.format(np.percentile(BW_KTrace,50)))
    print('第95分位数:{}'.format(np.percentile(BW_KTrace,95)))

    

    #ED  
    ED_Anpy=HQ_CR+r"\ED_A.npy"
    ED_Knpy=HQ_CR+r"\ED_K.npy"
    if os.path.exists(ED_Anpy)==False:
        ED_ATrace=attributeSampling.UniformSample(19,44)
        ED_KTrace=attributeSampling.UniformSample(5,6)
        print(ED_ATrace)
        print(ED_KTrace)
        np.save(ED_Anpy, ED_ATrace)
        np.save(ED_Knpy, ED_KTrace)
    ED_ATrace = np.load(ED_Anpy)
    ED_KTrace = np.load(ED_Knpy)
    # print('ED_A:')
    # print('第50分位数:{}'.format(np.percentile(ED_ATrace,50)))
    # print('第95分位数:{}'.format(np.percentile(ED_ATrace,95)))
    # print('ED_K:')
    # print('第50分位数:{}'.format(np.percentile(ED_KTrace,50)))
    # print('第95分位数:{}'.format(np.percentile(ED_KTrace,95)))


    
    #IR
    IRs_Anpy=HQ_CR+r"\IRs_A.npy"
    IRs_Knpy=HQ_CR+r"\IRs_K.npy"
    if os.path.exists(IRs_Anpy)==False:
        IRs_ATrace=attributeSampling.LognormalSample(50,76.5)
        IRs_KTrace=attributeSampling.LognormalSample(78,110.7)
        print(IRs_ATrace)
        print(IRs_KTrace)
        # print(IRd_ATrace)
        # print(IRd_KTrace)
        np.save(IRs_Anpy, IRs_ATrace)
        np.save(IRs_Knpy, IRs_KTrace)
    IRs_ATrace = np.load(IRs_Anpy)
    IRs_KTrace = np.load(IRs_Knpy)
    print('IRs_A:')
    print('第50分位数:{}'.format(np.percentile(IRs_ATrace,50)))
    print('第95分位数:{}'.format(np.percentile(IRs_ATrace,95)))
    print('IRs_K:')
    print('第50分位数:{}'.format(np.percentile(IRs_KTrace,50)))
    print('第95分位数:{}'.format(np.percentile(IRs_KTrace,95)))


    #SA
    SA_Anpy=HQ_CR+r"\SA_A.npy"
    SA_Knpy=HQ_CR+r"\SA_K.npy"
    if os.path.exists(SA_Anpy)==False:
        SA_KTrace=attributeSampling.SA_LognormalSample(8000,765.3)
        SA_ATrace=attributeSampling.SA_LognormalSample(16000,1530.6)
        print(SA_ATrace)
        print(SA_KTrace)
        np.save(SA_Anpy, SA_ATrace)
        np.save(SA_Knpy, SA_KTrace)
    SA_ATrace = np.load(SA_Anpy)
    SA_KTrace = np.load(SA_Knpy)
    print('SA_A:')
    print("负数：",SA_ATrace[SA_ATrace<0])
    print("负数：",SA_KTrace[SA_KTrace<0])
    print('第50分位数:{}'.format(np.percentile(SA_ATrace,50)))
    print('第95分位数:{}'.format(np.percentile(SA_ATrace,95)))
    print('SA_K:')
    print('第50分位数:{}'.format(np.percentile(SA_KTrace,50)))
    print('第95分位数:{}'.format(np.percentile(SA_KTrace,95)))


    #AF
    AF_Anpy=HQ_CR+r"\AF_A.npy"
    AF_Knpy=HQ_CR+r"\AF_K.npy"
    if os.path.exists(AF_Anpy)==False:
        AF_KTrace=attributeSampling.BetaSample(60,240)
        AF_ATrace=attributeSampling.BetaSample(70,930)
        print(AF_ATrace)
        print(AF_KTrace)
        np.save(AF_Anpy, AF_ATrace)
        np.save(AF_Knpy, AF_KTrace)
    AF_ATrace = np.load(AF_Anpy)
    AF_KTrace = np.load(AF_Knpy)
    print('AF_A:')
    print('第50分位数:{}'.format(np.percentile(AF_ATrace,50)))
    print('第95分位数:{}'.format(np.percentile(AF_ATrace,95)))
    print('AF_K:')
    print('第50分位数:{}'.format(np.percentile(AF_KTrace,50)))
    print('第95分位数:{}'.format(np.percentile(AF_KTrace,95)))

    #计算各健康风险指数
    attributeSampling.CalHI_CR(MC_SHP+r"\12PCDW_T",'12PCDW_T',BW_ATrace,BW_KTrace,ED_ATrace,ED_KTrace,IRs_ATrace,IRs_KTrace,SA_ATrace,SA_KTrace,AF_ATrace,AF_KTrace)
    #各风险指数出图
    show_cdf(HQ_CR+r"\HQ_CR_2.npy")

    #各影响因子栅格转GDM输入数据csv，输入至R语言代码运行GDM
    #步骤：1、创建一个渔网后ExtractMultiValuesToPoints.CopyFeature，，2、多值提取至点ExtractMultiValuesToPoints.ExtractMultiValuesToPoints，3、属性表导出为csv。arcpy、arcgis批量操作Export_ShpFieldValueToTxt.py，4、运行R代码







