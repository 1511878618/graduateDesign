import feature
protein = 'IVPNSVEQKHIQKED'
AAC= feature.CalculateAAComposition(protein)

DPC=feature.CalculateDipeptideComposition(protein)

CTD = feature.calculateCTD(protein)
CTDMatrix=[item for key,item in CTD.items()]
print(CTDMatrix)
#print(len(CTDMatrix))
#AAIndex = feature.CalculateAAIndex(protein)
#print('hi',AAIndex)
    
        
    
    