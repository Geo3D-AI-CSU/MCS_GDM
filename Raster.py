from genericpath import exists
from matplotlib.pyplot import step
from osgeo import ogr
import numpy
from osgeo import gdal, gdalconst
import sys
import math
import numpy as np
import numpy.ma as ma
import pymc as pm
import GetRasterEnt
import os
import attributeSampling as AS


np.set_printoptions(suppress=True)

def GetRasterValue(inputRaster,RasterName,RasterPoint):
    Data,geotransform,rows,cols,inputPrj,inputTrans =GetRasterEnt.GetRasterValue(inputRaster)
    print('proj',inputPrj)
    originX = geotransform[0]#影像左上角横坐标
    originY = geotransform[3]#
    print('originX',originX)
    print('originY', originY)

    pixelWidth = geotransform[1]#遥感图像的垂直空间分辨率或者南北方向上的像素分辨率
    pixelHeight = geotransform[5]#遥感图像的垂直空间分辨率
    centrolX = originX + pixelWidth / 2
    centrolY = originY + pixelHeight / 2

    print('Data:', Data)

    print('inshp:', RasterPoint)
    driver = ogr.GetDriverByName('ESRI Shapefile')

    ds_shp = driver.Open(RasterPoint)
    lyr = ds_shp.GetLayer(0)

    Feacount = lyr.GetFeatureCount()

    RasterValue_arr = [[None for j in range(cols)] for i in range(rows)]
    PointX_arr =  [[None for j in range(cols)] for i in range(rows)]
    PointY_arr =  [[None for j in range(cols)] for i in range(rows)]
    for i in range(Feacount):
        fea = lyr.GetFeature(i)
        feat_geom = fea.GetGeometryRef()
        geomX=feat_geom.GetX()
        geomY = feat_geom.GetY()
        #print("geomX:",geomX)
        #print("geomY:", geomY)
        field_index = fea.GetFieldIndex(RasterName)
        Value_j = fea.GetField(field_index)
        xOffset=int((geomX-centrolX)/pixelWidth)#列偏移数
        yOffset=int((geomY-centrolY)/pixelHeight)#行偏移数
        #print("xOffset:",xOffset)
        #print("yOffset:", yOffset)

        #栅格单元类别值
        RasterValue_arr[yOffset-1][xOffset-1]=Value_j

        #用于计算两单元距离
        PointX_arr[yOffset-1][xOffset-1] = geomX
        PointY_arr[yOffset-1][xOffset-1] = geomY


    print("len(RasterValue_arr)",len(RasterValue_arr))
    Raster_ValueNeighbor=[[None for j in range(cols)] for i in range(rows)]#每个像元点的邻域集合。将邻域数组按像元行列排列加入
    Raster_ArcDistanceNeighbor = [[None for j in range(cols)] for i in range(rows)]#对应的邻域像元与中心像元的反距离
    for i in range(len(PointX_arr)):
        for j in range(len(PointX_arr[i])):
            Value_neighbor = []
            Arcdistance_neighbor = []
            if (PointX_arr[i][j] is not None and PointX_arr[i][j]>0):
                target_X = PointX_arr[i][j]
                target_Y = PointY_arr[i][j]
                neighborStart_X = PointX_arr[i-2][j-2]
                neighborStart_Y = PointY_arr[i-2][j-2]
                #获取邻域起止行列号
                neighborStart_row=i-2
                neighborStart_col = j - 2

                for k in range(5):
                    neighbor_row = neighborStart_row + k
                    for l in range(5):
                        neighbor_col = neighborStart_col + l
                        #print("neighbor_row:",neighbor_row)
                        #print("neighbor_col:", neighbor_col)
                        #若超出了栅格边界，continue
                        if (neighbor_row <0 or neighbor_row>=rows or neighbor_col<0 or neighbor_col>=cols):
                            continue
                        elif (RasterValue_arr[neighbor_row][neighbor_col] is not None and RasterValue_arr[neighbor_row][neighbor_col] >0):
                            value_l = RasterValue_arr[neighbor_row][neighbor_col]
                            Value_neighbor.append(value_l)
                            neighbor_X= PointX_arr[neighbor_row][neighbor_col]
                            neighbor_Y =  PointX_arr[neighbor_row][neighbor_col]
                            #获取每个邻域像元与中心像元的反距离
                            ArcDistance_l=1/(math.sqrt((target_X - neighbor_X)**2+(target_Y - neighbor_Y)**2)+1)
                            Arcdistance_neighbor.append(ArcDistance_l)
            #print('len(Value_neighbor)',len(Value_neighbor))
            Raster_ValueNeighbor[i][j]=Value_neighbor
            Raster_ArcDistanceNeighbor[i][j]=Arcdistance_neighbor

    return rows,cols,RasterValue_arr,Raster_ValueNeighbor,Raster_ArcDistanceNeighbor,inputPrj,inputTrans

def getClassW(rows,cols,RasterValue_arr,Raster_ValueNeighbor,Raster_ArcDistanceNeighbor):
    Raster_WeightNeighbor=[[None for j in range(cols)] for i in range(rows)]#创建同样行列数的
    Raster_newClass = [[None for j in range(cols)] for i in range(rows)]
    for i in range(len(RasterValue_arr)):
        for j in range(len(RasterValue_arr[i])):
            W_neighbor=[]
            #计算反距离权重
            Arcdistance_totle=sum(Raster_ArcDistanceNeighbor[i][j])
            for k in range(len(Raster_ArcDistanceNeighbor[i][j])):
                if (Raster_ArcDistanceNeighbor[i][j][k] is not None):
                    W_j=Raster_ArcDistanceNeighbor[i][j][k]/Arcdistance_totle
                    #print("W_j:",W_j)
                    W_neighbor.append(W_j)
            Raster_WeightNeighbor[i][j]=W_neighbor
    for i in range(len(RasterValue_arr)):
        for j in range(len(RasterValue_arr[i])):
            #邻域内类别value去重
            NeighborClassArr = list(set(Raster_ValueNeighbor[i][j]))
            #创建字典类别：权重
            Class_Weight={x: {} for x in NeighborClassArr}

            for Class in NeighborClassArr:
                Weight=0
                for k in range(len(Raster_ValueNeighbor[i][j])):
                    if (Raster_ValueNeighbor[i][j] is not None):
                        if Raster_ValueNeighbor[i][j][k]==Class:
                            #print("Raster_WeightNeighbor[i][j][k]",Raster_WeightNeighbor[i][j][k])
                            Weight=round (Weight+Raster_WeightNeighbor[i][j][k],4)
                            Class_Weight[Class]=Weight
            if (len(NeighborClassArr)>1):
                print("ClassWeight:",Class_Weight)
            #模拟新栅格类别数据
            newClass=[None for n in range(100)]
            if (len(NeighborClassArr)==1):
                for n1 in range(100):
                    newClass[n1]=NeighborClassArr[0]
                    #print(newClass)
                Raster_newClass[i][j]=newClass#trace
                print("newclass:", newClass)
            elif (NeighborClassArr==[]):
                for n2 in range(100):
                    newClass[n2]=None
                    #print(newClass)
                Raster_newClass[i][j]=newClass  # trace
                print("newclass:",newClass)
            else:
                print("i{},j{}".format(i,j))
                trace_npy=r'E:\syc\myshengtai_soil\soilTest\Trace_npy\grid{}{}.npy'.format(i,j)

                newClass=BuildModel_U(Class_Weight)
                Raster_newClass[i][j]=newClass
    return Raster_newClass


def  BuildModel_U(ClassWeight):

    U=pm.Uniform('U', lower=0, upper=1)
    model = pm.Model([U])
    runner = pm.MCMC(model)
    # 进行10000次抽样，取后5000次
    runner.sample(iter=5000, burn=2000)
    UTrace = runner.trace("U")[:]
    #print("UTrace",UTrace)
    newTrace=[]
    for i in range(100):
        maxNum=0
        minNum=0
        #设置均匀分布，各类别所占比例区间
        for Class in ClassWeight.items():
            maxNum=maxNum+Class[1]
            UTrace=list(UTrace)
            if UTrace[i]>minNum and UTrace[i]<maxNum:
                #print("Class:",Class[0])
                Type=type(Class[0])
                UTrace[i]=Class[0]
                break
            minNum=maxNum
        newTrace.append(UTrace[i])
    return newTrace




def CreateRaster(arr,rows,cols,outPath,RasterName,inputPrj,inputTrans):
    Raster_arr=[[None for j in range(cols)] for i in range(rows)]
    for n in range(100):
        outRaster=outPath+"/"+RasterName+"_{}.tif".format(n+1)
        if not exists(outRaster):
            for i in range(rows):
                for j in range(cols):
                    Raster_arr[i][j]=arr[i][j][n]

            GetRasterEnt.arr2raster(Raster_arr, outRaster, inputPrj, inputTrans)


def RasterSampler_main(Raster_arr,RasterName_arr,Path):
    outdataPath = Path + r"\Raster_outdata5000"
    indataPath = Path + r"\Raster_indata5000"
    for i in range(len(Raster_arr)):
        outPath=outdataPath +'/'+RasterName_arr[i]
        inPath = indataPath + '/' + RasterName_arr[i]
        if os.path.exists(outPath):
            print("目录已存在。")
        else:
            os.makedirs(outPath)
            print("创建成功！")
        inputRaster=inPath+'/'+Raster_arr[i]

        RasterPoint=Path+'/'+r'fishNet_hunan5000\fishNet_hunan5000.shp'

        # print("RasterValue_arr:",RasterValue_arr)
        npy_Raster = outPath + '/' + RasterName_arr[i]+"5000_arr.npy"
        print('npy_Raster', npy_Raster)
        rows, cols, RasterValue_arr, Raster_ValueNeighbor, Raster_ArcDistanceNeighbor,inputPrj,inputTrans = GetRasterValue(inputRaster,RasterName_arr[i],
                                                                                                        RasterPoint)
        if not os.path.exists(npy_Raster):

            Raster_newClass = getClassW(rows, cols, RasterValue_arr, Raster_ValueNeighbor, Raster_ArcDistanceNeighbor)
            test_X = np.array(Raster_newClass)
            np.save(npy_Raster , test_X)
        RasterTrace = np.load(npy_Raster,allow_pickle=True)
        print("Raster_newClass:", RasterTrace)
        CreateRaster(RasterTrace, rows, cols, outPath, RasterName_arr[i],inputPrj,inputTrans)
















