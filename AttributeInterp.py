#!/usr/bin/env python
# _*_ coding:utf-8 _*_
__author__ = 'suyingcai'


import os
import sys
#print(sys.path)
import arcpy

from arcpy.sa import *
arcpy.env.overwriteOutput = True
# 设置输入点要素
#arcpy.env.extent =arcpy.Extent(37280398.2975199,2725633.3392,37823398.2975,3334633.3392)

#普通克里金插值
def interpolation(inshp,outshp,ditu,fieldName_arr):
    extent = arcpy.Describe(ditu).extent
    for i in range(0,100):
        #输入要素
        inFeatures =inshp+r"\Attribute12PCDW{}.shp".format(i+1)
        
        for j in range(len(fieldName_arr)):
            arcpy.env.extent =extent
            # Set Mask environment
            arcpy.env.mask = ditu
            outname=fieldName_arr[j]+"inter{}".format(i+1)
            # 设置输出结果
            outRaster = outshp+r'\{}interR'.format(fieldName_arr[j])+'\{}'.format(outname)
            
            # 设置输出字段
            field = fieldName_arr[j]
            # 设置输出栅格像元大小
            cellSize = 2160.83
            # 设置半变异函数属性
            kModelOrdinary = KrigingModelOrdinary("SPHERICAL", cellSize)
            # Check out the ArcGIS Spatial Analyst extension license 检查许可
            arcpy.CheckOutExtension("Spatial")
            # 执行克里金插值
            print('执行克里金插值，插值字段为：%s，输出像元大小为： %s'%(field,cellSize))
            # VARIABLE 12 是指 搜索半径设置，这里是arcgis默认的12
            outKriging = Kriging(inFeatures, field, kModelOrdinary, cellSize,"VARIABLE 12")
            print('保存插值结果至: %s' %outRaster)
            outKriging.save(outRaster)

#反距离权重插值
def IDWInter(inshp,outshp,ditu,fieldName_arr):
    
    for i in range(100):
        # Set local variables
        inPointFeatures =inshp+r"\12PCDW_T{}.shp".format(i+1)
        print('inPointFeatures:',inPointFeatures)
        for fieldName in fieldName_arr:
            print("fieldname:",fieldName)

            #Set Mask environment
            arcpy.env.mask = ditu
            outname=fieldName +"{}".format(i+1)
            # 设置输出结果
            outRaster = outshp+r'\{}'.format(fieldName)+'\{}.tif'.format(outname)
            print('outRaster',outRaster)
            
            zField = fieldName
            cellSize =1000
            power = 2
            #searchRadius = RadiusVariable(15, 199326)

            # Check out the ArcGIS Spatial Analyst extension license
            #arcpy.CheckOutExtension("Spatial")

            # Execute IDW
            #outIDW = Idw(inPointFeatures, zField, cellSize, power, searchRadius)
            #print('保存插值结果至: %s' %outRaster)
            # Save the output 
            #outIDW.save(outRaster)

            # Check out the ArcGIS 3D Analyst extension license
            arcpy.CheckOutExtension("3D")
             
            # Execute IDW
            arcpy.Idw_3d(inPointFeatures, zField, outRaster, cellSize,power)
            print("插值完成")
            
        
            
            
        




    

def CalMeanSTD(interPath,outpath,fieldName_arr):
    for n in range(len(fieldName_arr)):
        inshp_arr=[]
        inPath=interPath+r"\{}_PR".format(fieldName_arr[n])
        outPath= outpath+r'\{}_PR'.format(fieldName_arr[n])
        
        for i in range(100):
            FeatureName=fieldName_arr[n]+"{}M".format(i+1)
            inFeature = inPath+'\{}'.format(FeatureName)
            print("inFeature:",inFeature)
            if not os.path.exists(inPath):
                print("路径不存在")
            inshp_arr.append(inFeature)
        # Execute CellStatistics
        outSTDname=outPath+r"\{}STD".format(fieldName_arr[n])
        outMEANname=outPath+r"\{}MEAN".format(fieldName_arr[n])
        print("outSTDname:",outSTDname)
        print("outMEANname:",outMEANname)
        # Check out the ArcGIS Spatial Analyst extension license
        arcpy.CheckOutExtension("Spatial")
        outShpSTD = CellStatistics(inshp_arr, "STD","NODATA")
        outShpMEAN = CellStatistics(inshp_arr, "MEAN","NODATA")
        # Save the output
        
        outShpSTD.save(outSTDname)
        outShpMEAN.save(outMEANname)
        


def SetProj(inPath,outPath,ditu,fieldName_arr):
    #extent = arcpy.Describe(ditu).extent
    #arcpy.env.extent =extent
    #Set Mask environment
    #arcpy.env.mask = ditu
    for i in range(100):
        # Set local variables
        for fieldName in fieldName_arr:
            print("fieldname:",fieldName)
            rasterName=fieldName +"IDW{}".format(i+1)
            # 设置输出结果
            inRaster = inPath+r'\{}_PR'.format(fieldName)+'\{}'.format(rasterName)
            print('inRaster',inRaster)
            try:
                # set local variables
                in_dataset = inRaster#"forest.shp"
                
                # get the coordinate system by describing a feature class
                dsc = arcpy.Describe(ditu)
                sr = arcpy.SpatialReference(4546)#CGCS2000 111E
                #coord_sys = dsc.spatialReference
                
                # run the tool
                arcpy.DefineProjection_management(in_dataset, sr)
                
                # print messages when the tool runs successfully
                print(arcpy.GetMessages(0))
                
            except arcpy.ExecuteError:
                print(arcpy.GetMessages(2))
                
            except Exception as ex:
                print(ex.args[0])
            #ExtractMask(inRaster,ditu,outPath,fieldName,i)
        
        
            
            
    
    
    
def ExtractMask(RasterPath,Shp,OutPath,fieldName_arr):
    for i in range(100):
        # Set local variables
        for fieldName in fieldName_arr:
            rasterName=fieldName +"IDW{}".format(i+1)
            # 设置输入栅格
            inRaster = RasterPath+r'\{}_PR'.format(fieldName)+'\{}'.format(rasterName)
            print('inRaster',inRaster)
            # Set local variables
            #inRaster = Raster
            inMaskData = Shp
            OutName=fieldName+"{}M".format(i+1)
            OutRaster = OutPath+r'\{}_PR'.format(fieldName)+'\{}'.format(OutName)
            print("OutRaster",OutRaster)

            # Check out the ArcGIS Spatial Analyst extension license
            arcpy.CheckOutExtension("Spatial")

            # Execute ExtractByMask
            outExtractByMask = ExtractByMask(inRaster, inMaskData)

            # Save the output 
            outExtractByMask.save(OutRaster)
            print(u"保存掩膜结果成功")
            
        
    


if __name__ == "__main__":
    Path=r'D:\syc\lessons\master\meeting\progress\data\test'
    inshp=Path+r'\testAttribute'
    outPath=r""
    ditu=Path+r"\行政边界_省级.shp"
    fieldName_arr=["NO2","PM2_5"]

    #shengtai
    FieldName_arr=['Nemerow']
    #FieldName_arr=['PCDW_Hg']
    dataPath=r"E:\syc\myshengtai_soil"
    indataPath=dataPath+r"\indata"
    outdataPath=dataPath+r"\soilTest"
    PCDWAttr=outdataPath+r'\12PCDWAttrSampling'
    PCDWPoint=outdataPath+r'\12PCDWSampling'
    PCDWPoint_10=outdataPath+r'\12PCDW_T'
    hunan=r'E:\syc\myshengtai_soil\soilTest\12PCDW_2000.shp'


    #普通克里金插值
    #interpolation(PCDWAttr,outdataPath,hunan,FieldName_arr)

    #IDW插值
    IDWInter(PCDWPoint_10,outdataPath,hunan,FieldName_arr)
    #SetProj(outdataPath,outdataPath,hunan,FieldName_arr)
    #ExtractMask(outdataPath,hunan,outdataPath,FieldName_arr)
    #CalMeanSTD(outdataPath,outdataPath,FieldName_arr)
    
