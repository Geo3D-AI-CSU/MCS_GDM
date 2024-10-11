import numpy
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
import sys
import os

#from Raster import RasterName_arr
os.environ['PROJ_LIB'] = r'D:\anaconda3\envs\pymc\Library\share\proj'
os.environ['GDAL_DATA'] = r'D:\anaconda3\envs\pymc\Library\share'

import threading
import numpy as np

def GetRasterValue(inputRaster):
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    ds = gdal.Open(inputRaster, 1)
    if ds is None:
        print('Could not open ' + inputRaster)
        sys.exit(1)
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    bands = ds.RasterCount
    print("cols",cols)
    print("rows", rows)
    print("bands", bands)

    geotransform = ds.GetGeoTransform()#仿射矩阵
    inputPrj = ds.GetProjection()
    inputTrans = ds.GetGeoTransform()
    print(inputPrj)
    print(inputTrans)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(-999)
    print(band.GetNoDataValue())
    
    data= band.ReadAsArray(0, 0, cols-1, rows-1)
    print(data)

    return data,geotransform,rows,cols,inputPrj,inputTrans


def GetRasterValue100(startCount,endCount,RasterName,RasterPath,Dataset):
    for i in range(startCount,endCount,1):
        Raster = RasterPath + '/'+RasterName+"_{}.tif".format(i + 1)
        #print(Raster)
        Data,geotransform,rows,cols,inputprj,inputtrans =GetRasterValue(Raster)
        Dataset.append(Data)



def arr2raster(arr, raster_file, prj=None, trans=None):
    """
    将数组转成栅格文件写入硬盘
    :param arr: 输入的mask数组 ReadAsArray()
    :param raster_file: 输出的栅格文件路径
    :param prj: gdal读取的投影信息 GetProjection()，默认为空
    :param trans: gdal读取的几何信息 GetGeoTransform()，默认为空
    :return:
    """
    arr = numpy.array(arr)
    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(raster_file, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)

    if prj:
        dst_ds.SetProjection(prj)
    if trans:
        dst_ds.SetGeoTransform(trans)

    # 将数组的各通道写入图片
    outband=dst_ds.GetRasterBand(1)
    outband.WriteArray(arr)
    outband.SetNoDataValue(-999)

    dst_ds.FlushCache()
    outband.GetStatistics(0, 1)
    dst_ds = None
    print("successfully convert array to raster")


#转换DataSet组织方式
def TransformDataSet(Dataset):
    DataSet = []
    for rows in range(len(Dataset[0])):#行数/FeaCount
        grid_row=[]
        for cols in range(len(Dataset[0][0])):#列数/FieldCount
            grid_col = []#该像元的100个取值/shp数1-100
            for i in range(len(Dataset)):
                GridValue=Dataset[i][rows][cols]
                print('GridValues:', GridValue)
                grid_col.append(GridValue)
            grid_row.append(grid_col)
        DataSet.append(grid_row)
    print('TransformDataSet:',DataSet)
    return DataSet


# 计算分类信息熵
def calcShannonEnt(dataset):
    # 初始化变量
    SumEnt = 0.0  # 所有特征的总熵
    Ent_arr = []  # 存储每个特征熵值的列表

    dataSet = np.array(dataset)  # 将数据集转换为NumPy数组
    print(len(dataSet))  # 打印数据集的行数

    # 遍历数据集中的每一行
    for row in range(len(dataSet)):
        Ent_row = []  # 存储每个行中每个值的熵的列表

        # 遍历行中每个特征的值
        for values in dataSet[row]:
            gridEnt = 0.0  # 初始化当前特征的熵

            # 创建特征中唯一值的集合
            value_list = set([values[i] for i in range(values.shape[0])])
            num_of_values = len(values)

            print(values)
            print("value_list:", value_list)

            # 计算特征中每个唯一值的熵
            for value in value_list:
                p = float(values[values == value].shape[0]) / values.shape[0]  # 计算值的概率
                logp = np.log2(p)  # 计算以2为底的对数
                gridEnt -= p * logp  # 更新特征的熵

            print(gridEnt)
            Ent_row.append(gridEnt)  # 将特征的熵追加到行的熵列表中

        Ent_arr.append(Ent_row)  # 将行的熵列表追加到整体熵列表中

    return Ent_arr


def EntMethod(Path,RasterNameArr):
    for RasterName in RasterNameArr:
        RasterPath=Path+'/' + RasterName
        outEntPath=RasterPath

        inputPrj = 'PROJCS["CGCS2000_3_Degree_GK_Zone_37",GEOGCS["GCS_China_Geodetic_Coordinate_System_2000",\
        DATUM["China_2000",SPHEROID["CGCS2000",6378137,298.257222101,AUTHORITY["EPSG","1024"]],AUTHORITY["EPSG","1043"]],PRIMEM["Greenwich",0],\
        UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",111],\
        PARAMETER["scale_factor",1],PARAMETER["false_easting",37500000],\
        PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
        #inputTrans =(37275398.29751994, 5000.0, 0.0, 3339633.339203578, 0.0, -5000.0)#可用，但范围不对
        inputTrans =(37275398.29751994, 5000.0, 0.0, 3339633.339203578, 0.0, -5000.0)

        npy_Raster = RasterPath+ '/' + RasterName + "Ent_arr.npy"
        outEntRaster=outEntPath+ '/' + RasterName + "Ent.tif"
        print('npy_Raster', npy_Raster)
        if not os.path.exists(npy_Raster):
            processes = []
            Dataset = []
            dataset1 = []
            dataset2 = []
            dataset3 = []
            dataset4 = []


            driver = gdal.GetDriverByName('HFA')
            driver.Register()

            #多线程读取数据后将返回的数组拼接
            poc1 = threading.Thread(target=GetRasterValue100,\
                                    args=(0, 25, RasterName, RasterPath,  dataset1))
            poc2 = threading.Thread(target=GetRasterValue100,\
                                    args=(25, 50, RasterName, RasterPath,  dataset2))
            poc3 = threading.Thread(target=GetRasterValue100,\
                                    args=(50, 75, RasterName, RasterPath,  dataset3))
            poc4 = threading.Thread(target=GetRasterValue100,\
                                    args=(75, 100, RasterName, RasterPath, dataset4))
            for p in (poc1, poc2, poc3, poc4):
                processes.append(p)
                p.start()
            for p in processes:
                p.join()
            Dataset = dataset1 + dataset2 + dataset3 + dataset4
            print('sumDataset:', Dataset)
            print('RowCount:', len(Dataset[0]))
            print('ColCount:', len(Dataset[0][0]))
            print('RasterCount:', len(Dataset))
            #转换数组维度
            DataSet = TransformDataSet(Dataset)
            #计算信息熵
            Ent_arr = calcShannonEnt(DataSet)
            test_X = np.array(Ent_arr)
            np.save(npy_Raster, test_X)
        RasterTrace = np.load(npy_Raster, allow_pickle=True)
        #生成信息熵栅格
        arr2raster(RasterTrace, outEntRaster, prj=inputPrj, trans=inputTrans)  #
        # print("PointEnt:",PointEnt_arr[0])
        # print("PointEnt_arr:",PointEnt_arr)
        # print("SumEnt:",SumEnt)





