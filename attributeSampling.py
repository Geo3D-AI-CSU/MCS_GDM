from multiprocessing import Value
from symbol import parameters
from osgeo import ogr, osr
#import matplotlib.pyplot as plt
# %matplotlib inline
import numpy as np
import os
# import matplotlib.pyplot as plt
import pandas as pd
import pymc as pm
from osgeo import ogr,gdal
import math
import main as M
import GetRasterEnt



#点插值属性的不确定性。
def GetFeas(inshp,fieldName_arr,fieldsType_arr):
    print('inshp:', inshp)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")

    gdal.SetConfigOption('SHAPE_ENCODING', "CP936")  # 属性表支持中文字段
    ds = driver.Open(inshp)
    lyr = ds.GetLayer()
    prjfile = open(inshp.replace('.shp', '.prj'), 'w')
    Feacount = lyr.GetFeatureCount()
    geoX_arr=[]
    geoY_arr=[]
    geoFieldsValue_arr=[]

    lyrDefn = lyr.GetLayerDefn()
    feature0=lyr.GetFeature(0)
    featureDefn = feature0.GetDefnRef()
    if fieldsType_arr == []:
        for fieldName in fieldName_arr:
            field_index = lyr.FindFieldIndex(fieldName, True)
            fieldDefn = lyrDefn.GetFieldDefn(field_index)
            field_type = fieldDefn.GetTypeName()
            if field_type=='Real':
                field_type=ogr.OFTReal
            elif field_type=='String':
                field_type=ogr.OFTString
            print(field_type)
            fieldsType_arr.append(field_type)

    for i in range(Feacount):
        fea = lyr.GetFeature(i)
        feat_geom = fea.GetGeometryRef()
        geoX=feat_geom.GetX()
        geoY=feat_geom.GetY()
        geoX_arr.append(geoX)
        geoY_arr.append(geoY)
        FieldsValue_arr = []
        for fieldName in fieldName_arr:
            field_index = fea.GetFieldIndex(fieldName)
            Value = fea.GetField(field_index)
            FieldsValue_arr.append(Value)


        geoFieldsValue_arr.append(FieldsValue_arr)

    return geoX_arr,geoY_arr,geoFieldsValue_arr,fieldsType_arr

                      
def BuildModel(AttributesArr,sigma_arr):
    
    for i in range(len(AttributesArr)):#每个geom
        AttributeTraceArr = []
        print("geoAttrArr:", AttributesArr[i])

        for j in range(len(AttributesArr[i])):#每个字段
            if (AttributesArr[i][j]==0):               
                AttributesArr[i][j]=sigma_arr[j] #避免抽样出现负值，若初始值为0，则将初始值 更改为sigma大小
            print("geo{},Field{}".format(i,j))
            trace_npy=M.MC_ATTRIBUTE+r'\attribute_Test{}{}.npy'.format(i,j)
            if os.path.exists(trace_npy)==False:
                m=AttributesArr[i][j]
                v = sigma_arr[j]
                # mu=math.log((m**2)/(math.sqrt(v+m**2)))
                # sigma=math.sqrt(math.log(v/(m**2)+1))
                mu = np.log(m / np.sqrt(1 + (v / m)**2))                
                sigma = np.sqrt(np.log(1 + (v / m)**2))
                tau=1/(sigma*sigma)
                Attribute = pm.Normal("attribute{}{}".format(i,j), mu=mu, tau=tau)
                model = pm.Model([Attribute])
                runner = pm.MCMC(model)
                # 进行5000次抽样，取后2000次
                runner.sample(iter=5000, burn=2000)
                AttributeTrace = runner.trace("attribute{}{}".format(i,j))[:]
                test = np.array(AttributeTrace)
                np.save(trace_npy, test)
    geoAttributesTrace_arr=[]
    for i0 in range(len(AttributesArr)):#每个geom
        AttributeTraceArr = []
        #print("geoAttrArr:", AttributesArr[i])
        for j0 in range(len(AttributesArr[i])):#每个字段   
            trace_npy=M.MC_ATTRIBUTE+r'\attribute_Test{}{}.npy'.format(i0,j0) 
            trace = np.load(trace_npy)
            AttributeTraceArr.append(trace)
        geoAttributesTrace_arr.append(AttributeTraceArr)  
    return geoAttributesTrace_arr


def GB15618(k,PH,waterField):
    if k == 0:  # Cr 铬
        if PH <= 5.5:
            if waterField == -1:  # -1代表为水田
                background = 250
            else:
                background = 150
        elif PH <= 6.5:
            if waterField == -1:  
                background = 250
            else:
                background = 150
        elif PH <= 7.5:
            if waterField == -1:  
                background = 300
            else:
                background = 200
        else:
            if waterField == -1:  
                background = 350
            else:
                background = 250
    if k == 1:  # Cd 镉
        if PH <= 5.5:
            if waterField == -1:  
                background = 0.3
            else:
                background = 0.3
        elif PH <= 6.5:
            if waterField == -1:  
                background = 0.4
            else:
                background = 0.3
        elif PH <= 7.5:
            if waterField == -1:  
                background = 0.6
            else:
                background = 0.3
        else:
            if waterField == -1:  
                background = 0.8
            else:
                background = 0.6
    if k == 2:  # As 砷
        if PH <= 5.5:
            if waterField == -1:  
                background = 30
            else:
                background = 40
        elif PH <= 6.5:
            if waterField == -1:  
                background = 30
            else:
                background = 40
        elif PH <= 7.5:
            if waterField == -1:  
                background = 25
            else:
                background = 30
        else:
            if waterField == -1:  
                background = 20
            else:
                background = 25
    if k == 3:  # Pb 铅
        if PH <= 5.5:
            if waterField == -1:  
                background = 80
            else:
                background = 70
        elif PH <= 6.5:
            if waterField == -1:  
                background = 100
            else:
                background = 90
        elif PH <= 7.5:
            if waterField == -1:  
                background = 140
            else:
                background = 120
        else:
            if waterField == -1:  
                background = 240
            else:
                background = 170
    if k == 4:  # Hg 汞
        if PH <= 5.5:
            if waterField == -1:  
                background = 0.5
            else:
                background = 1.3
        elif PH <= 6.5:
            if waterField == -1:  
                background = 0.5
            else:
                background = 1.8
        elif PH <= 7.5:
            if waterField == -1:  
                background = 0.6
            else:
                background = 2.4
        else:
            if waterField == -1:  
                background = 1.0
            else:
                background = 3.4
    return background



#各元素非致癌风险
def HQ(kind,CS,BW_A,BW_K,ED_A,ED_K,IRs_A,IRs_K,SA_A,SA_K,AF_A,AF_K):
    #暴露途径：口摄 s，皮肤接触 d；人群：A表示为成人，K表示为儿童；
    # IRs_A=100#摄入率
    # IRs_K=200#
    #AFd=0.2
    TF=10**-6
    EF_A=350#
    EF_K=350
    # ED_A=24
    # ED_K=6
    # BW_A=63
    # BW_K=29
    AT_A=365*ED_A
    AT_K=365*ED_K
    AT_C=365*70#致癌
    # SA_A=1.6*10**4
    # SA_K=2800
    if kind=='Cr':
        ABS=0.001
        RfDs=3*10**-3
        RfDd=6*10**-5
    elif kind=='Cd':
        ABS=0.001
        RfDs=1*10**-3
        RfDd=1*10**-5
    elif kind=='As':
        ABS=0.03
        RfDs=3*10**-4
        RfDd=1.23*10**-4
    elif kind=='Pb':
        ABS=0.001
        RfDs=3.5*10**-3
        RfDd=5.25*10**-4
    elif kind=='Hg':
        ABS=0.001
        RfDs=3*10**-4
        RfDd=2.1*10**-5

    #非致癌
    ADDs_A=CS*IRs_A*TF*EF_A*ED_A/(BW_A*AT_A)
    ADDs_K=CS*IRs_K*TF*EF_K*ED_K/(BW_K*AT_K)
    ADDd_A=CS*AF_A*TF*SA_A*ABS*EF_A*ED_A/(BW_A*AT_A)
    ADDd_K=CS*AF_K*TF*SA_K*ABS*EF_K*ED_A/(BW_K*AT_K)
    #致癌
    ADDs_AC=CS*IRs_A*TF*EF_A*ED_A/(BW_A*AT_C)
    ADDs_KC=CS*IRs_K*TF*EF_K*ED_K/(BW_K*AT_C)
    ADDd_AC=CS*AF_A*TF*SA_A*ABS*EF_A*ED_A/(BW_A*AT_C)
    ADDd_KC=CS*AF_K*TF*SA_K*ABS*EF_K*ED_A/(BW_K*AT_C)
    

    HQs_A=ADDs_A/RfDs
    HQs_K=ADDs_K/RfDs
    HQd_A=ADDd_A/RfDd
    HQd_K=ADDd_K/RfDd

    return ADDs_A,ADDs_K,ADDd_A,ADDd_K,ADDs_AC,ADDs_KC,ADDd_AC,ADDd_KC,HQs_A,HQs_K,HQd_A,HQd_K

#致癌
def CR(kind,ADDs_AC,ADDs_KC,ADDd_AC,ADDd_KC):
    if kind=='Cr':
        SFs=5*10**-1
        SFd=20  ##
    elif kind=='Cd':
        SFs=5.01*10**-1
        SFd=20
    elif kind=='As':
        SFs=1.5
        SFd=3.66
    elif kind=='Pb':#无用，后续不讨论Pb的致癌风险，故在后续代码中将其统一赋值为-999
        SFs=1
        SFd=1
    elif kind=='Hg':#无用，后续不讨论Pb的致癌风险，故在后续代码中将其统一赋值为-999
        SFs=1
        SFd=1
    CRs_A=ADDs_AC*SFs
    CRs_K=ADDs_KC*SFs
    CRd_A=ADDd_AC*SFd
    CRd_K=ADDd_KC*SFd
    return CRs_A,CRs_K,CRd_A,CRd_K

def CalHI_CR(shpPath,shpName,BW_ATrace,BW_KTrace,ED_ATrace,ED_KTrace,IRs_ATrace,IRs_KTrace,SA_ATrace,SA_KTrace,AF_ATrace,AF_KTrace):
    FieldNameArr=['Cr','Cd','As','Pb','Hg']    #元素字段  
    fieldsType_arr=[]
    FieldsValue_arr=[]
    #读取100个shp各要素各字段的值
    for i0 in range(100):
        inshp = shpPath + '/' + shpName + r"{}.shp".format(i0 + 1)
        #print('inshp:', inshp)
        geoX_arr,geoY_arr,geoFieldsValue_arr,fieldsType_arr=GetFeas(inshp, FieldNameArr,fieldsType_arr)
        FieldsValue_arr.append(geoFieldsValue_arr)
    Value100Arr=[]   
    #每个shp 
    for i1 in range(100):
        # 计算每个点        
        ValueArr=[]
        for j in range(len(FieldsValue_arr[i1])):
            FeaIndexValueArr=[]
            # 计算每种重金属的HQ、CR
            for k1 in range(len(FieldNameArr)):
                ADDs_A,ADDs_K,ADDd_A,ADDd_K,ADDs_AC,ADDs_KC,ADDd_AC,ADDd_KC,HQs_A,HQs_K,HQd_A,HQd_K=HQ(FieldNameArr[k1],FieldsValue_arr[i1][j][k1],BW_ATrace[i1],BW_KTrace[i1],ED_ATrace[i1],ED_KTrace[i1],IRs_ATrace[i1],IRs_KTrace[i1],SA_ATrace[i1],SA_KTrace[i1],AF_ATrace[i1],AF_KTrace[i1])
                HI_A=HQs_A+HQd_A
                HI_K=HQs_K+HQd_K
                if FieldNameArr[k1] =='Pb' or FieldNameArr[k1] =='Hg' :
                    CR_A=-999
                    CR_K=-999
                else: 
                    CRs_A,CRs_K,CRd_A,CRd_K=CR(FieldNameArr[k1],ADDs_AC,ADDs_KC,ADDd_AC,ADDd_KC)
                    CR_A=CRs_A+CRd_A
                    CR_K=CRs_K+CRd_K
                FeaIndexValueArr.append(HQs_A)
                FeaIndexValueArr.append(HQs_K)
                FeaIndexValueArr.append(HQd_A)
                FeaIndexValueArr.append(HQd_K)
                FeaIndexValueArr.append(CRs_A)
                FeaIndexValueArr.append(CRs_K)
                FeaIndexValueArr.append(CRd_A)
                FeaIndexValueArr.append(CRd_K)
                FeaIndexValueArr.append(HI_A)
                FeaIndexValueArr.append(HI_K)
                FeaIndexValueArr.append(CR_A)
                FeaIndexValueArr.append(CR_K)               
            ValueArr.append(FeaIndexValueArr)
        Value100Arr.append(ValueArr)
    npy_Value100Arr=M.HQ_CR+r"\Value100Arr_2.npy"
    # print("Value100arr shape:",Value100Arr)
    if os.path.exists(npy_Value100Arr)==False:
        np.save(npy_Value100Arr, Value100Arr,allow_pickle=True)
    Value100Arr = np.load(npy_Value100Arr)
    npy_value=M.HQ_CR+r"\HQ_CR_2.npy"
    if os.path.exists(npy_value)==False:
        TransformData = Value100Arr.transpose(1, 2, 0)
        print("TransformData shape:",TransformData.shape)
        np.save(npy_value, TransformData,allow_pickle=True)
    
    TransformData = np.load(npy_value)
    print(TransformData.shape)
    return




def calNemerow(shpPath,shpName,Si):
    FieldNameArr=['Cr','Cd','As','Pb','Hg','PH',u"水田"]
    if Si == 'GB15816':
        fieldName_arr=['CrIndex','CdIndex','AsIndex','PbIndex','HgIndex','Nemerow']
    fieldsType_arr=[]
    FieldsValue_arr=[]
    for i0 in range(100):
        inshp = shpPath + '/{}/'.format(shpName) + shpName + r"{}.shp".format(i0 + 1)
        #inshp = shpPath + '/soilTest/' + shpName + r".shp"
        #print('inshp:', inshp)
        geoX_arr,geoY_arr,geoFieldsValue_arr,fieldsType_arr=GetFeas(inshp, FieldNameArr,fieldsType_arr)
        FieldsValue_arr.append(geoFieldsValue_arr)
        
    IndexValue100Arr=[]
    for i1 in range(100):
        # 计算每个点
        ValueArr=[]
        for j in range(len(FieldsValue_arr[i1])):
            FeaIndexValueArr=[]
            # 计算每种重金属的单因素污染指数
            PH = FieldsValue_arr[i1][j][5]
            waterField = FieldsValue_arr[i1][j][6]

            for k in range(5):
                if Si=='GB15816':
                    background=GB15618(k, PH, waterField)
                value=FieldsValue_arr[i1][j][k]
                SingleIndex=value/background
                FeaIndexValueArr.append(SingleIndex)
            maxSingle=np.max(FeaIndexValueArr)
            meanSingle=np.sum(FeaIndexValueArr)/len(FeaIndexValueArr)
            NemerowIndex=np.sqrt((maxSingle*maxSingle+meanSingle*meanSingle)/2)
            FeaIndexValueArr.append(NemerowIndex)
            ValueArr.append(FeaIndexValueArr)
        IndexValue100Arr.append(ValueArr)
    Trans_value=M.MC_SHP+'/'+shpName+r'\Trans_value.npy'
    if os.path.exists(Trans_value)==False:
        TansformData=GetRasterEnt.TransformDataSet(IndexValue100Arr)
        #print("geoAttributesTrace_arr",geoAttributesTrace_arr)
        test = np.array(TansformData)
        np.save(Trans_value, test)

    TansformData = np.load(Trans_value)
    WriteField(shpPath, shpName, TansformData, fieldName_arr,fieldsType_arr[0:6])



#给图层中多个字段同时赋值，适用于土壤点抽样后数据赋值
def WriteField(shpPath,shpName,geoFieldsValue_arr,fieldName_arr,field_typeArr):
    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")

    gdal.SetConfigOption('SHAPE_ENCODING', "CP936")  # 属性表支持中文字段
    driver = ogr.GetDriverByName('ESRI Shapefile')

    for n in range(0,100):
        #inshp = shpPath +'/'+shpName+ r"{}.shp".format(n + 1)
        inshp =  shpPath + '/{}/'.format(shpName) + shpName + r"{}.shp".format(n + 1)
        print('inshp:', inshp)
        
        ds = driver.Open(inshp,1)
        lyr = ds.GetLayer()

        lyrDefinition = lyr.GetLayerDefn()
        fieldCount = lyrDefinition.GetFieldCount()
        for j in range(len(fieldName_arr)):
            # print(j)          
            print("fieldName:", fieldName_arr[j])
            fieldIndex = lyrDefinition.GetFieldIndex(fieldName_arr[j])
            if fieldIndex < 0:
                Field = ogr.FieldDefn(fieldName_arr[j], field_typeArr[j])
                lyr.CreateField(Field, 1)
                print("创建字段完成！")
                # print("value:",geoFieldsValue_arr[i][j][n])
            ndim = np.array(geoFieldsValue_arr).ndim
            # print("数组维度：", ndim)
            for i in range(len(geoFieldsValue_arr)):
                print(i)
                feature=lyr.GetFeature(i)
                if ndim >2:
                    #土壤的对数按正态分布抽样后，反对数变换
                    value=math.exp(geoFieldsValue_arr[i][j][n])                      
                    feature.SetField(fieldName_arr[j],round(value,4))
                    lyr.SetFeature(feature)
                    
                else:
                    feature.SetField(fieldName_arr[j], geoFieldsValue_arr[i][j])
                    lyr.SetFeature(feature)
           



# 反距离权重插值,也可使用arcpy、arcgis软件批量处理
def IDWInter(inshp, outshp, fieldName_arr, outName_arr):
    for j in range(len(fieldName_arr)):

        # 设置输出结果
        outPath = outshp + r'\{}'.format(outName_arr[j])

        if os.path.exists(outPath):
            print("目录已存在。")
        else:
            os.makedirs(outPath)
            print("创建成功！")
        for i in range(100):
            # Set local variables
            inPointFeatures = inshp + r"\12PCDW_T{}.shp".format(i + 1)
            print('inPointFeatures:', inPointFeatures)
            outname = fieldName_arr[j] + "{}.tif".format(i + 1)
            outRaster = outPath + '\{}'.format(outname)
            print('outRaster', outRaster)

            zField = fieldName_arr[j]
            print("zField",zField)
            """
            	idw空间插值
            	:param output_file:插值结果
            	:param point_station_file: 矢量站点数据
                spatFilter=(minX,minY,maxX,maxY),spatFilter=(37280398.2975199,2725633.33920358,37823398.2975199,37823398.2975199),
                outputBounds=(37280398.2975199,2725633.3392,37823398.2975,3334633.3392)
                spatFilter=(37280399.1274,2725201.5327,37825399.1274,3335201.5327)
            	:return:
            	"""
            # 代码调用
        
            opts = gdal.GridOptions(algorithm='invdist:power=2:smoothing=0.2:radius1=0.0:radius2=0.0:angle=0.0:max_points=0:min_points=0:nodata=0.0',\
                                    format="GTiff", outputType=gdal.GDT_Float32, width=543,height=609,outputSRS='EPSG:4525',\
                                    outputBounds=(37280398.2975199,2725633.3392,37823398.2975,3334633.3392),zfield=zField,noData=-9999)
            gdal.Grid(destName=outRaster, srcDS=inPointFeatures,options=opts)
            print("插值完成")
            
        
def ExtractMask(RasterPath,Shp,OutPath,fieldName_arr):
    for i in range(0,100):
        # Set local variables
        for fieldName in fieldName_arr:
            #rasterName=fieldName +"{}.tif".format(i+1)
            rasterName=fieldName +"Ent.tif"#.format(i+1)
            # 设置输入栅格
            inRaster = RasterPath+r'\{}'.format(fieldName)+'\{}'.format(rasterName)
            print('inRaster',inRaster)
            # Set local variables
            #inRaster = Raster
            inMaskData = Shp
            OutName=fieldName+"{}_m.tif".format(i+1)
            OutRaster = OutPath+r'\{}'.format(fieldName)+'\{}'.format(OutName)
            print("OutRaster",OutRaster)

            # Execute ExtractByMask
            ds = gdal.Warp(OutRaster , inRaster, format='GTiff',\
                           cutlineDSName=inMaskData, dstNodata=-9999, dstSRS='EPSG:4546',srcSRS='EPSG:4525')
            print(u"保存掩膜结果成功")
        
            
        


def main(inshp,outshp,fieldName_arr,sigma_arr):
    FieldName_arr = ['Cr', 'Pb', 'Cd', 'As', 'Hg','PH']
    #FieldName_arr = [u'水田']
    fieldsTypeArr=[]
    #获取原始字段取值
    geoX_arr,geoY_arr,geoFieldsValue_arr,fieldsType_arr=GetFeas(inshp,fieldName_arr,fieldsTypeArr)
    npy_value=outshp+r"\traceValue0.2.npy"
    if os.path.exists(npy_value)==False:

        geoAttributesTrace_arr=BuildModel(geoFieldsValue_arr,sigma_arr)

        #print("geoAttributesTrace_arr",geoAttributesTrace_arr)
        test = np.array(geoAttributesTrace_arr)
        np.save(npy_value, test)
    trace = np.load(npy_value)
    WriteField(M.MC_SHP,'12PCDW_T',trace, FieldName_arr, fieldsType_arr)
    


# AF Beta抽样
def BetaSample(alpha,beta):
    
    B=pm.Beta('B', alpha=alpha,beta=beta)#alpha+beta越大则抽样数值越集中，alpha/(alpha+beta)越小越靠近0
    model = pm.Model([B])
    runner = pm.MCMC(model)
    # 进行5000次抽样，取后3000次
    runner.sample(iter=5000, burn=2000)
    BTrace = runner.trace("B")[:]   
    return BTrace

#
def LognormalSample(mu,sigma):
    
    L=pm.Lognormal('L', mu=mu,tau=1/(sigma*sigma))#
    model = pm.Model([L])
    runner = pm.MCMC(model)
    # 进行5000次抽样，取后3000次
    runner.sample(iter=5000, burn=2000)
    LTrace = runner.trace("L")[:]
    print("小于1:",LTrace[LTrace<1])
    if len(LTrace[LTrace<1])>0:
        LTrace[LTrace<1]=np.random.choice(LTrace[LTrace>1],len(LTrace[LTrace<1]))#随机选择进行替换
    print("小于1_new:",LTrace[LTrace<1])
    print(len(LTrace))
    LTrace=[math.log(x) for x in LTrace]
    print('第50分位数:{}'.format(np.percentile(LTrace,50)))
    print('第95分位数:{}'.format(np.percentile(LTrace,95)))
    return LTrace

def SA_LognormalSample(m,v):
    # mu=math.log((m**2)/(math.sqrt(v+m**2)))
    # sigma=math.sqrt(math.log(v/(m**2)+1))
    mu = np.log(m / np.sqrt(1 + (v / m)**2))
    sigma = np.sqrt(np.log(1 + (v / m)**2))
    tau=1/(sigma*sigma)
    L=pm.Normal('L', mu=mu,tau=tau)#
    model = pm.Model([L])
    runner = pm.MCMC(model)
    # 进行5000次抽样，取后3000次
    runner.sample(iter=5000, burn=2000)
    LTrace = runner.trace("L")[:]  
    LTrace=[math.exp(x) for x in LTrace]
    print('第50分位数:{}'.format(np.percentile(LTrace,50)))
    print('第95分位数:{}'.format(np.percentile(LTrace,95)))
    return LTrace

def UniformSample(low,up):
    U=pm.Uniform('U', lower=low, upper=up)#
    model = pm.Model([U])
    runner = pm.MCMC(model)
    # 进行5000次抽样，取后3000次
    runner.sample(iter=5000, burn=2000)
    UTrace = runner.trace("U")[:] 
    print(UTrace)
    return UTrace


    
