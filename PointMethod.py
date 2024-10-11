

# %matplotlib inline
import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
import pymc as pm
from osgeo import ogr
import os


def GetPointFea(inshp,shpName,outPath):
    print('inshp:', inshp)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(inshp)
    lyr = ds.GetLayer()

    Feacount = lyr.GetFeatureCount()
    print("Feacount:", Feacount)
    i = 0
    geomX_arr = []
    geomY_arr = []
    for i in range(Feacount):
        fea = lyr.GetFeature(i)
        feat_geom = fea.GetGeometryRef()
        geomX = feat_geom.GetX()
        geomY = feat_geom.GetY()
        print("geomX:", geomX)
        print("geomY:", geomY)
        geomX_arr.append(geomX)
        geomY_arr.append(geomY)
        print("####################################################")
    sigma_x = 10
    sigma_y = 10
    npy_X = outPath + '/' + shpName + '_X2.npy'
    npy_Y = outPath + '/' + shpName + '_Y.npy'
    if not os.path.exists(npy_X):
        PointXTrace_arr, PointYTrace_arr, runner = build_model(geomX_arr, geomY_arr, sigma_x, sigma_y)
        test_X = np.array(PointXTrace_arr)
        np.save(outPath + '/' + shpName + '_X', test_X)
    trace_X = np.load(npy_X)
    if not os.path.exists(npy_Y):
        test_Y = np.array(PointYTrace_arr)
        np.save(outPath + '/' + shpName + '_Y', test_Y)
    trace_Y = np.load(npy_Y)

    return trace_X,trace_Y

def CreatPointShp(FeaXTrace_arr, FeaYTrace_arr, Type,shpname,outpath):
    for i in range(100):  # 
        outname = shpname + "{}.shp".format(i + 1)
        outshp = outpath + "/" + outname
        print("outshp:", outshp)
        driver = ogr.GetDriverByName('ESRI Shapefile')
        data_source = driver.CreateDataSource(outshp)
        #srs = osr.SpatialReference()
        #srs.ImportFromEPSG(4326)  # Xian80)
        layer = data_source.CreateLayer("ceshi", None, Type)

        createShp(layer, i, FeaXTrace_arr, FeaYTrace_arr, Type)
        data_source.FlushCache()


def Point_Method(PointShp_arr,shpName_Arr, Path):
    Type = ogr.wkbPoint
    for i in range(len(PointShp_arr)):
        outPath=Path+'/'+shpName_Arr[i]
        if os.path.exists(outPath):
            print("已存在")
        else:
            os.makedirs(outPath)
            print("创建成功")

        trace_X,trace_Y=GetPointFea(PointShp_arr[i],shpName_Arr[i], outPath)
        CreatPointShp(trace_X,trace_Y, Type, shpName_Arr[i], outPath)





def createShp(layer, i, XTrace, YTrace, Type):
    print("进入模型！")
    print(XTrace)
    if (Type == ogr.wkbPoint):

        feature = ogr.Feature(layer.GetLayerDefn())

        print("获取到feature")
        for j in range(len(XTrace)):
            print(j)
            point = ogr.Geometry(Type)
            point.AddPoint(XTrace[j][i], YTrace[j][i])

            feature.SetGeometry(point)

            layer.CreateFeature(feature)
            print("创建 完成")

        feature.Destroy()




def build_model(XList,YList,sigma_x,sigma_y):
    params=[]
    for i in range(len(XList)):#建立所选点参数
        X=pm.Normal('x{}'.format(i),mu=XList[i],tau=1/sigma_x**2)
        Y=pm.Normal('y{}'.format(i),mu=YList[i],tau=1/sigma_y**2)
        params.append(X)
        params.append(Y)
    print("params:",len(params))
    model = pm.Model(params)
    runner = pm.MCMC(model)
    #进行10000次抽样，取后5000次
    runner.sample(iter =5000, burn=2000)
    xTrace_arr,yTrace_arr=[],[]
    for j in range(len(XList)):
        xTrace=runner.trace("x{}".format(j))[:]
        yTrace=runner.trace("y{}".format(j))[:]
        xTrace_arr.append(xTrace)
        yTrace_arr.append(yTrace)

    return xTrace_arr,yTrace_arr,runner
