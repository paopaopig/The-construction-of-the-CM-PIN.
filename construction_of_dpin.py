import numpy as np
import pandas as pd
species = 'human'
dataset = 'HDIP'

SPIN = np.load(r"model_detection/%s/%s/spin/spin.npy"%(species,dataset))
nodes_names = np.load(r"model_detection/%s/%s/nodes_name.npy"%(species,dataset))
        
nodes_expression_value = np.load(r"model_detection/%s/%s/gene_value_profile.npy"%(species,dataset))



def get_threshold(k,after_combine_matrix):
    sigma = np.std(after_combine_matrix,axis = 1,ddof=1)
    mean = np.mean(after_combine_matrix,axis = 1)
    th = mean+(k*sigma)
    return th
k = 0.5


yuzhi = get_threshold(k, nodes_expression_value)


for i in range(len(nodes_names)):
    for j in range(nodes_expression_value.shape[1]):
        if yuzhi[i] < nodes_expression_value[i,j]:
            nodes_expression_value[i,j] = 1
        else:
            nodes_expression_value[i,j] = 0
            

DPIN = np.zeros((len(nodes_names), len(nodes_names)))

for i in range(len(nodes_names)):
    for j in range(i+1,len(nodes_names)):
        if (nodes_expression_value[i,:]*nodes_expression_value[j,:]).sum()>=1:
            DPIN[i,j] = DPIN[j,i] = 1  
        
DPIN = DPIN * SPIN
print(sum(sum(DPIN))/2)
np.save(r'model_detection/%s/%s/dpin/dpin-k=0.5.npy'%(species,dataset),DPIN)     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    