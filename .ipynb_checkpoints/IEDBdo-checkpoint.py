import pandas as pd
import copy 

from Bio import SeqIO
from Bio import * 

import numpy as np 

"""
这个包只用于处理IEDB下载的数据，简略处理，方便调用。
"""
AALteter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]


# 筛选数据
##筛选出非bces
def deduplicateFromTotalData(smallData,TotalData):
    """
    input:
            smallData:是一个小数据集 
            TotalData:是一个大数据，里面包含里小数据集和另外的一些数据
    output:
            outData:是大数据集中除了小数据集以外的其他数据集
    """
    index = [i for i in TotalData.index if i not in smallData.index]
    outData = TotalData.loc[index]
    return outData

##去除包含非20种基本氨基酸的多肽序列
def deleteWrongAmino(Data):
    '''
    input:
        Data：是一个df，其中Description是多肽序列，去掉包含有非基本氨基酸的多肽序列
        AALteter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    outout: 
    '''
    data = copy.deepcopy(Data)
    #筛选是否存在氨基酸不是20种基本氨基酸中的一类
    c= []
    for index,values in data['Description'].items():
        for amino in values:
            if amino not in AALteter:
                c.append(index)
                print('IEDBID为{}，的多肽序列包含有未知氨基酸{}'.format(index,amino))
                break 
    if c:
        data = data.drop(labels=c,axis=0)
    else:
        pass 
    return data

## 去除修饰
def delModifications(Data):
    """
    input:data是一个df，其中Description是多肽序列，且修饰是 AAALPG + S12 这样的格式
    outpu:去除了所有的修饰，保留蛋白质序列
    """
    data =copy.deepcopy(Data)
    data['Description'] = data['Description'].apply(lambda x: x.split(' + ')[0])
    return data 
#保存数据
##保存数据为fasta格式
def saveSeqToFasta(df,filename):
    data = copy.deepcopy(df)
    seqdata = data['Description'].apply(lambda x : Seq.Seq(x))
    seqRecords = [SeqRecord.SeqRecord(seq,id=str(ID),description='') for ID,seq in seqdata.items()]
    SeqIO.write(seqRecords,filename,'fasta')

##读取fasta格式的数据
def readSeqFromeFasta(filename):
    seqData = np.array([str(seq.seq) for seq in list(SeqIO.parse(filename,'fasta'))]).astype('object')
    return seqData

#数据处理pepline

def eatIEDBData(positiveData,totalData,dataSetName,returnSeq=True):
    """
    input:
        positiveData：是IEDB上勾选实验验证为阳性的数据集
        totalData：是IEDB上在同样的筛选条件下不勾选实验验证为阳性后得到的包含阳性和非阳性的数据集
        datasetName：是处理后生成的两个fasta文件的名称，该名称由：dataSetName+(non)BCEs.fastq 组成
        returnSeq：是否返回两个数据集，若返回则是np.array类型的两个数组，dtype为object
    output:
        在本路径下保存处理后生成的fasta文件
        returnSeq：返回两个数据集，若返回则是np.array类型的两个数组，dtype为objectrn
    调用示例：
         Total = pd.read_csv('人食物过敏.csv',header=1,index_col=0)
         BCEs = pd.read_csv('人食物过敏相关BCEs.csv',header=1,index_col=0)
         a,b = eatIEDBData(BCEs,Total,'人食物过敏数据集')
         
    
    """
    #得到非bces
    nonBCEs = deduplicateFromTotalData(positiveData,totalData)
    #去除修饰
    BCEs = delModifications(positiveData)
    nonBCEs = delModifications(nonBCEs)
    #删除包含非基本氨基酸的多肽序列
    BCEs = deleteWrongAmino(BCEs)
    nonBCEs = deleteWrongAmino(nonBCEs)
    #保存数据
    BCEsFileName = dataSetName+'_'+'BCEs.fasta'
    saveSeqToFasta(BCEs,BCEsFileName)
    nonBCEsFileName = dataSetName+'_'+'nonBCEs.fasta'
    saveSeqToFasta(nonBCEs,nonBCEsFileName)
    #读取数据
    if returnSeq:
        return readSeqFromeFasta(BCEsFileName),readSeqFromeFasta(nonBCEsFileName)