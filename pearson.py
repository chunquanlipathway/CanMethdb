# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 19:23:25 2018

@author: Administrator
"""
import os,sys
from scipy import stats
import pandas as pd
import numpy as np
import datetime

def txt2dataframe(files):
    '''
    功能：读取txt文件，存到dataframe数据类型中
    files：txt文件，读取为dataframe类型数据，并返回
    '''
    f = pd.read_table(files,engine='python',encoding='Utf-8',index_col=False,sep="\t")
    #打印内存使用量
    print '------------------------内存使用量------------------------'
    print files+' : '
    print f.info(memory_usage='deep')
    print '--------------------------------------------------------'
    return f

def getvalue(name,namelist,namedf):
    '''
    功能：获取查询条目在dataframe中的行，并返回，查询不到返回空
    name:查询名称
    namelist:要查询名称的列表
    namedf:存储列表及相应数值的dataframe
    '''
    if name in namelist:
        return namedf.ix[namelist.index(name)].values.tolist()
    else:
        return []
def writetxt(of,context):
    with open(of,'a+') as f:
        for i in context:
            f.write(str(i)+'\n')
            
def main(matrix_f,meth_f,pair_f,output,matrix_title,meth_title):
    '''
    功能：求皮尔逊相关系数
    matrix_f:基因表达谱
    meth_f:甲基化普
    pair_f:关系对
    matrix_title:表达谱第一列列名
    meth_title:甲基化谱第一列列名
    '''
    #读取txt文件为frame类型
    matrix_df = txt2dataframe(matrix_f)
    meth_df = txt2dataframe(meth_f)
    #去除表达谱中，‘。’后的字符串，存储到列表
    gene_list = [gene.split('.')[0] for gene in matrix_df[matrix_title].values.tolist()]
    #获取甲基化普位点列表
    methlocaloc_list = meth_df[meth_title].values.tolist()
    #以追加写的模式，打开输出文件
    header = meth_title+'\t'+matrix_title+'\t'+'r-value'+'\t'+'p-value'+'\n'
    of = open(output,'a+')
    of.write(header)
    #逐行读取关系对文件
    with open(pair_f) as f:
        #获取表头
        header = f.readline()
        #记录行数
        number = 1 
        #循环关系对文件的每一行
        for i in f:
            #打印行数
            if number % 10000 == 0:
                print str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+' : '+str(number)
            number = number + 1
            #对每行文件去除换行符，并以\t进行分割
            line = i.replace('\r','').strip('\n').split('\t')
            #基因名
            methloca = line[0]
            #甲基化位点
            genesymb = line[1]
            #从dataframe中获取基因相对应的每个样本表达值
            x = getvalue(genesymb,gene_list,matrix_df)[1:]
            #从dataframe中获取甲基化位点对应的值
            y = getvalue(methloca,methlocaloc_list,meth_df)[1:]
            #如果x,y有一个为空，则跳出循环，不求取皮尔逊相关系数
            if len(x)==0 or len(y)==0:
                continue
            else:
                #print len(x),len(y)
                #计算皮尔逊相关系数，返回值，第一位是相关系数，第二位是p值
                #print x,y
                #
                #print len(x),len(y)
                #writetxt('gene.txt',x)
                #writetxt('meth.txt',y)
                cor = stats.pearsonr(x[1:],y[1:])
                print cor
                of.write(methloca+'\t'+genesymb+'\t'+str(cor[0])+'\t'+str(cor[1])+'\n')
                of.flush()
                #break
    of.close()

if __name__ == '__main__':
    if len(sys.argv)<6:
        print "参数数量不对!"
        exit()
    matrix_f = sys.argv[1] #'pearson/RNA_comm_no0.txt'
    meth_f = sys.argv[2] #'pearson/Methy_comm_nona.txt'
    pair_f = sys.argv[3] #'pearson/all_pair.bed'
    output = sys.argv[4] #'pearson/result2.txt'
    #判断输入文件是否存在，存在则继续执行，不存在则结束程序
    for i in [matrix_f,meth_f,pair_f]:
        if os.path.exists(i):
            pass
        else:
            print i+" ==> No such file!"
            exit()
    #如果输出文件存在，则结束程序。
    if os.path.exists(output):
        print output+' ==> File exists!'
        exit()
    else:
        matrix_title = sys.argv[5] #"V1"
        meth_title = sys.argv[6] #"Composite.Element.REF"
        main(matrix_f,meth_f,pair_f,output,matrix_title,meth_title)
