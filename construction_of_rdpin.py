import numpy as np
import pandas as pd

species = 'human'
dataset = 'HDIP'

DPIN = np.load(r'model_detection/%s/%s/dpin/dpin-k=0.5.npy'%(species,dataset))   
nodes_names = np.load(r"model_detection/%s/%s/nodes_name.npy"%(species,dataset))
  

sublocation = pd.read_csv(r'model_detection/%s/%s/subloc.csv'%(species,dataset))
gene_name=np.array(sublocation.iloc[:,0])
subcell_nuc=np.array(sublocation.iloc[:,1:])# 11个亚细胞室

k=0
nodes_nucleus = []
for i in range(len(nodes_names)):
    flag = 0
    if nodes_names[i] in gene_name:
        ind = np.argwhere(gene_name == nodes_names[i])
        flag += 1
        k=k+1
        nodes_nucleus.append(subcell_nuc[ind[0,0],:])
    if flag == 0:
         nodes_nucleus.append(np.zeros(11))

nodes_nucleus = np.array(nodes_nucleus)


RDPIN = np.zeros((len(nodes_names), len(nodes_names)))
for i in range(len(nodes_names)):
    for j in range(i+1, len(nodes_names)):
        if (nodes_nucleus[i,:]*nodes_nucleus[j,:]).sum()>=1:
            RDPIN[i,j] = RDPIN[j,i] = 1

RDPIN = RDPIN*DPIN
print(sum(sum(RDPIN))/2)
np.save(r'model_detection/%s/%s/rdpin/rdpin-k=0.5.npy'%(species,dataset),RDPIN)     
            