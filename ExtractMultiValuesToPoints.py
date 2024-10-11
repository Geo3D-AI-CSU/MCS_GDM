#!/usr/bin/env python
# _*_ coding:utf-8 _*_
__author__ = 'suyingcai'


import arcpy
from arcpy import *
from arcpy.sa import *
from arcpy import env
from glob import glob            #方便遍历shp文件
import sys,os


env.workspace = path



def ExtractMultiValuesToPoints(Raster_outdata,OutNameArr)
    for i in range(0,100):
        inRasterlist=[]
        for j in range(len(OutNameArr)):
            raster=Raster_outdata+"/"+Raster_name[j]+"/{}.tif".format(Raster_name[j])
            print(raster)
            RasterList=[raster,OutNameArr[j]]
            inRasterlist.append(RasterList)
        shp=path+r'\12PCDW_T{}.shp'.format(i+1)
        print(inRasterlist)
        ExtractMultiValuesToPoints(shp,inRasterlist, "NONE")
        
    print('ok')
    


def CopyFeature(temFea,outPath):
    
    for i in range(0,100):
        outFea=outPath+r"\fishNet_hunan{}".format(i+1)
        outLayer=arcpy.CopyFeatures_management(temFea,outFea)
    print("done!")

    
if __name__ == "__main__":

    fishNet_hunan5000=r'E:\syc\soil\fishNet_hunan5000\fishNet_hunan5000.shp'
    fishNet_hunan5000Path=r'E:\syc\soil\fishNet_hunan5000'
    path=r"E:\syc\myshengtai_soil\soilTest\12PCDW_T"
    Raster_outdata=path+"/Raster_indata5000"
    OutNameArr=['Slope','Aspect','ST','WSL','pH','GDP','PD','MA','LULC','DRW','DR','DW','TEM','PRE','SOC','dzt','Nemerow']
    #arcmap创建一个渔网后，copy
    CopyFeature(fishNet_hunan5000,fishNet_hunan5000Path)
    #多值提取至点
    ExtractMultiValuesToPoints(Raster_outdata,OutNameArr)
