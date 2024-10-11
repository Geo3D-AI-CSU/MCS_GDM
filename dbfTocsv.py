
#!/usr/bin/env python
# _*_ coding:utf-8 _*_
__author__ = 'suyingcai'


import os
import sys
#print(sys.path)
import arcpy

def Export_ShpFieldValueToTxt(src_data, dst_data, field_list):
    arcpy.env.workspace = r'C:/Users/Administrator/Documents/ArcGIS/Default.gdb'
    fp = open(dst_data, 'w')
    fp.write('FID,Slope,Aspect,ST,WSL,pH,MA,LULC,PRE,SOC,GDP,PD,DRW,DR,DW,TEM,Lith,Nemerow')
    #fp.write('Nemerow,X1')
    fp.write('\n')
    with arcpy.da.SearchCursor(src_data, field_list) as cursor:
        for row in cursor:
            # print(row[0])
            if (row[0] == ''):
                continue
            for j in range(len(field_list)):
                if field_list[j]=='Nemerow':
                    fp.write(str(row[j]))
                #print(str(row[j]))
                else:
                    fp.write(str(row[j])+",")
            fp.write('\n')
    fp.close()

field_list=['FID','Slope','Aspect','ST','WSL','pH','MA','LULC','PRE','SOC','GDP','PD','DRW','DR','DW','TEM','dzt','Nemerow']
#field_list=['Nemerow','X1']
for i in range(100):
    print(i)
    shp=r'E:\syc\soil\fishNet_hunan5000\fishNet_hunan{}.shp'.format(i+1).decode('gbk')
    print(shp)
    csv=r'E:\syc\soil\fishNet_NemerowCSV5000\fishNet_NemerowCSV{}.csv'.format(i+1).decode('gbk')
    Export_ShpFieldValueToTxt(shp, csv, field_list)
#Export_ShpFieldValueToTxt(r'E:\syc\soil\fishNet_hunan5000\fishNet_hunanYuan.shp', r'E:\syc\soil\fishNet_NemerowCSV5000\fishNet_NemerowYuanCSV.csv', field_list)    
print('ok')
