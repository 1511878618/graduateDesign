


"""
可以计算多肽序列的特征值！
    Amino Acid Composition AAC
    Dipeptide Composition DPC
    Chain-transition-distribution CTD


"""

AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
#############################################################################################
def CalculateAAComposition(ProteinSequence):

    """
    ########################################################################
    Calculate the composition of Amino acids

    for a given protein sequence.

    Usage:

    result=CalculateAAComposition(protein)

    Input: protein is a pure protein sequence.

    Output: result is a dict form containing the composition of

    20 amino acids.
    ########################################################################
    """
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    LengthSequence=len(ProteinSequence)
    Result={}
    for i in AALetter:
        Result[i]=round(float(ProteinSequence.count(i))/LengthSequence,3)
    return Result


def CalculateDipeptideComposition(ProteinSequence):
    LenthSequence = len(ProteinSequence)-1
    Result={}
    for i in AALetter:
        for j in AALetter:
            Dipeptide=i+j
            Result[Dipeptide]=round(float(ProteinSequence.count(Dipeptide))/LenthSequence,3)
    return Result




AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

Hydrophobicity={'1':'RKEDQN','2':'GASTPHY','3':'CLVIMFW'} 
#'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

NormalizedVDWV={'1':'GASTPD','2':'NVEQIL','3':'MHKFRYW'}
#'1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

Polarity={'1':'LIFWCMVY','2':'CPNVEQIL','3':'KMHFRYW'}
#'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

Charge={'1':'KR','2':'ANCQGHILMFPSTWYV','3':'DE'}
#'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'}
#'1'stand for Helix; '2'stand for Strand, '3' stand for coil

SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND','3':'MPSTHY'}
#'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

Polarizability={'1':'GASDT','2':'CPNVEQIL','3':'KMHFRYW'}
#'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)


##You can continuely add other properties of AADs to compute descriptors of protein sequence.

AATProperty=(Hydrophobicity,NormalizedVDWV,Polarity,Charge,SecondaryStr,SolventAccessibility,Polarizability)

AATPropertyName=['Hydrophobicity','NormalizedVDWV','Polarity','Charge','SecondaryStr','SolventAccessibility','Polarizability']
AATProperty = dict([(AATPropertyName[i],AATProperty[i]) for i in range(len(AATPropertyName))])

##################################################################################################




def stringtoNum(proteinSequence,AAPName):
    """
    ###############################################################################################
    Tranform the protein sequence into the string form such as 32123223132121123.

    Usage:

    result=StringtoNum(protein,AAProperty)

    Input: protein is a pure protein sequence.

    AAPName is a str , amino acids property such as Polarizability.

    Output: result is a string such as 123321222132111123222
    ###############################################################################################
    """
    for index,aa in AATProperty[AAPName].items():
        for a in aa:
            try:
                proteinSequence=proteinSequence.replace(a,index)
            except:
                pass
    return proteinSequence

def calculateComposition(proteinSequence,AAPName):
    """
    ###############################################################################################
    A method used for computing composition descriptors.

    Usage:

    result=CalculateComposition(protein,AAProperty,AAPName)

    Input: protein is a pure protein sequence.

    AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

    AAPName is a string used for indicating a AAP name.

    Output: Result是一个字典，包含各个class的含量百分比
    ###############################################################################################
    """
    hardProteinSequence = stringtoNum(proteinSequence,AAPName)#转化成数字
    Result={}
    Lenth=len(hardProteinSequence)
    #计算1,2,3的C值
    Result[AAPName+'C'+'1']=round(float(hardProteinSequence.count('1'))/Lenth,3)
    Result[AAPName+'C'+'2']=round(float(hardProteinSequence.count('2'))/Lenth,3)
    Result[AAPName+'C'+'3']=round(float(hardProteinSequence.count('3'))/Lenth,3)
    return Result

def calculateTransition(proteinSequence,AAPName):
    """
    input: 
        proteinSequence:是一个蛋白质序列的字符串
        AAPName:是要计算的氨基酸物理化学性质，是一个字符串
    output:
        Result：是包含了计算基于输入的物化性质的 transition 百分比的字典
    """
    hardProteinSequence = stringtoNum(proteinSequence,AAPName)#转化成数字
    Result={}
    Lenth=len(hardProteinSequence)-1 
    Result[AAPName+'T'+'12']=round(float(hardProteinSequence.count('12')+hardProteinSequence.count('21'))/Lenth,3)
    Result[AAPName+'T'+'13']=round(float(hardProteinSequence.count('13')+hardProteinSequence.count('31'))/Lenth,3)
    Result[AAPName+'T'+'23']=round(float(hardProteinSequence.count('23')+hardProteinSequence.count('32'))/Lenth,3)
    return Result

def calculateDistribution(proteinSequence,AAPName):
    """
    input: 
        proteinSequence:是一个蛋白质序列的字符串
        AAPName:是要计算的氨基酸物理化学性质，是一个字符串
    output:
        Result：是包含了计算基于输入的物化性质的 distribution 百分比的字典
    """
    hardProteinSequence = stringtoNum(proteinSequence,AAPName)#转化成数字
    Result={}
    Lenth=len(hardProteinSequence)

    for i in ['1','2','3']:
        DList=[]
        num = hardProteinSequence.count(i)#属性i的氨基酸个数
        #选出属性为i的氨基酸，并用元组标上每一个氨基酸在原序列的位置n+1
        for n in range(Lenth):
            Class=hardProteinSequence[n]
            if Class == i:
                DList.append((n+1,i))
        try:
            Result[AAPName+'D'+i+'001']=round(float(DList[0][0])/Lenth,3)
            Result[AAPName+'D'+i+'025']=round(DList[int(float(num*0.25))-1][0]/Lenth,3)
            Result[AAPName+'D'+i+'050']=round(DList[int(float(num*0.5))-1][0]/Lenth,3)
            Result[AAPName+'D'+i+'075']=round(DList[int(float(num*0.75))-1][0]/Lenth,3)
            Result[AAPName+'D'+i+'100']=round(float(DList[-1][0])/Lenth,3)
        except: 
            Result[AAPName+'D'+i+'001']=0.0
            Result[AAPName+'D'+i+'025']=0.0
            Result[AAPName+'D'+i+'050']=0.0
            Result[AAPName+'D'+i+'075']=0.0
            Result[AAPName+'D'+i+'100']=0.0
        
    return Result
    
def CalculateCTD(proteinSequence,AAPNameList=AATPropertyName):
    """
    input:
        AAPNameList is a list of ['Hydrophobicity','NormalizedVDWV',
                                'Polarity','Charge','SecondaryStr',
                                'SolventAccessibility','Polarizability']
        proteinSequence is a pure protein sequence str
    
    output:
        result is a dict with CTD 
    """
    result = {}
    if isinstance(AAPNameList,list):
        for AAPName in AAPNameList:
            result.update(calculateComposition(proteinSequence,AAPName))
            result.update(calculateTransition(proteinSequence,AAPName))
            result.update(calculateDistribution(proteinSequence,AAPName))

        return result
    else:
        print("请输入包含氨基酸物化性质的一个列表，如['Hydrophobicity'] or ['Hydrophobicity','NormalizedVDWV']")
        



AAIdex_8features={'BLAM930101':{'A':0.96,'R':0.77,'N':0.39,'D':0.42,'C':0.42,'Q':0.8,'E':0.53,'G':0.0,'H':0.57,'I':0.84,'L':0.92,'K':0.73,'M':0.86,'F':0.59,'P':-2.5,'S':0.53,'T':0.54,'W':0.58,'Y':0.72,'V':0.63},'BIOV880101':{'A':16.0,'R':-70.0,'N':-74.0,'D':-78.0,'C':168.0,'Q':-73.0,'E':-106.0,'G':-13.0,'H':50.0,'I':151.0,'L':16.0,'K':-70.0,'M':-74.0,'F':-78.0,'P':168.0,'S':-73.0,'T':-106.0,'W':-13.0,'Y':50.0,'V':151.0},'MAXF760101':{'A':1.43,'R':1.18,'N':0.64,'D':0.92,'C':0.94,'Q':1.22,'E':1.67,'G':0.46,'H':0.98,'I':1.04,'L':1.36,'K':1.27,'M':1.53,'F':1.19,'P':0.49,'S':0.7,'T':0.78,'W':1.01,'Y':0.69,'V':0.98},'TSAJ990101':{'A':89.3,'R':190.3,'N':122.4,'D':114.4,'C':102.5,'Q':146.9,'E':138.8,'G':63.8,'H':157.5,'I':163.0,'L':163.1,'K':165.1,'M':165.8,'F':190.8,'P':121.6,'S':94.2,'T':119.6,'W':226.4,'Y':194.6,'V':138.2},'NAKH920108':{'A':9.36,'R':0.27,'N':2.31,'D':0.94,'C':2.56,'Q':1.14,'E':0.94,'G':6.17,'H':0.47,'I':13.73,'L':9.36,'K':0.27,'M':2.31,'F':0.94,'P':2.56,'S':1.14,'T':0.94,'W':6.17,'Y':0.47,'V':13.73},'CEDJ970104':{'A':7.9,'R':4.9,'N':4.0,'D':5.5,'C':1.9,'Q':4.4,'E':7.1,'G':7.1,'H':2.1,'I':5.2,'L':8.6,'K':6.7,'M':2.4,'F':3.9,'P':5.3,'S':6.6,'T':5.3,'W':1.2,'Y':3.1,'V':6.8},'LIFS790101':{'A':0.92,'R':0.93,'N':0.6,'D':0.48,'C':1.16,'Q':0.95,'E':0.61,'G':0.61,'H':0.93,'I':1.81,'L':-0.37,'K':0.33,'M':-0.3,'F':-0.38,'P':0.19,'S':0.12,'T':0.03,'W':-0.33,'Y':-0.29,'V':-0.29}}

""" 
8个特征分别是：
BLAM930101, BIOV880101, MAXF760101, TSAJ990101, NAKH920108, CEDJ970104, LIFS790101, and MIYS990104

"""

def CalculateAAIndex(proteinSequence):
    """
    Properitys= {'BLAM930101': {'A': 0.96, ' R': 0.77, ' N': 0.39, ' D': 0.42, ' C': 0.42, ' Q': 0.8, ' E': 0.53, ' G': 0.0, ' H': 0.57, ' I': 0.84, ' L': 0.92, ' K': 0.73, ' M': 0.86, ' F': 0.59, ' P': -2.5, ' S': 0.53, ' T': 0.54, ' W': 0.58, ' Y': 0.72, ' V': 0.63}, 'BIOV880101': {'A': 16.0, ' R': -70.0, ' N': -74.0, ' D': -78.0, ' C': 168.0, ' Q': -73.0, ' E': -106.0, ' G': -13.0, ' H': 50.0, ' I': 151.0, ' L': 16.0, ' K': -70.0, ' M': -74.0, ' F': -78.0, ' P': 168.0, ' S': -73.0, ' T': -106.0, ' W': -13.0, ' Y': 50.0, ' V': 151.0}, 'MAXF760101': {'A': 1.43, ' R': 1.18, ' N': 0.64, ' D': 0.92, ' C': 0.94, ' Q': 1.22, ' E': 1.67, ' G': 0.46, ' H': 0.98, ' I': 1.04, ' L': 1.36, ' K': 1.27, ' M': 1.53, ' F': 1.19, ' P': 0.49, ' S': 0.7, ' T': 0.78, ' W': 1.01, ' Y': 0.69, ' V': 0.98}, 'TSAJ990101': {'A': 89.3, ' R': 190.3, ' N': 122.4, ' D': 114.4, ' C': 102.5, ' Q': 146.9, ' E': 138.8, ' G': 63.8, ' H': 157.5, ' I': 163.0, ' L': 163.1, ' K': 165.1, ' M': 165.8, ' F': 190.8, ' P': 121.6, ' S': 94.2, ' T': 119.6, ' W': 226.4, ' Y': 194.6, ' V': 138.2}, 'NAKH920108': {'A': 9.36, ' R': 0.27, ' N': 2.31, ' D': 0.94, ' C': 2.56, ' Q': 1.14, ' E': 0.94, ' G': 6.17, ' H': 0.47, ' I': 13.73, ' L': 9.36, ' K': 0.27, ' M': 2.31, ' F': 0.94, ' P': 2.56, ' S': 1.14, ' T': 0.94, ' W': 6.17, ' Y': 0.47, ' V': 13.73}, 'CEDJ970104': {'A': 7.9, ' R': 4.9, ' N': 4.0, ' D': 5.5, ' C': 1.9, ' Q': 4.4, ' E': 7.1, ' G': 7.1, ' H': 2.1, ' I': 5.2, ' L': 8.6, ' K': 6.7, ' M': 2.4, ' F': 3.9, ' P': 5.3, ' S': 6.6, ' T': 5.3, ' W': 1.2, ' Y': 3.1, ' V': 6.8}, 'LIFS790101': {'A': 0.92, ' R': 0.93, ' N': 0.6, ' D': 0.48, ' C': 1.16, ' Q': 0.95, ' E': 0.61, ' G': 0.61, ' H': 0.93, ' I': 1.81, ' L': -0.37, ' K': 0.33, ' M': -0.3, ' F': -0.38, ' P': 0.19, ' S': 0.12, ' T': 0.03, ' W': -0.33, ' Y': -0.29, ' V': -0.29}}
    input:
            proteinsequence: 是一个蛋白质序列， str
    output:
            result: 是一个包含20*8=160 个键值对的 字典 包含每一个氨基酸理化性质的值的均值
    
    """
    #AAIdex_8features
    result={}

    for ProperityName,Properity in AAIdex_8features.items():
        sequence = proteinSequence
        for aminoAcid,aminoAcidvalue in Properity.items():
            result[ProperityName+aminoAcid] = round(float(sequence.count(aminoAcid)*aminoAcidvalue)/len(sequence),3)
    return result



    