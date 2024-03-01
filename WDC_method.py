import numpy as np
import pandas as pd

def wdc(Dis, nodes_name, gene_value):
    gene_value = gene_value.astype("float64")
    degree = np.sum(Dis,axis = 1)
    
    
    Ecc = np.zeros((len(nodes_name),len(nodes_name)))
    P = np.zeros((len(nodes_name),len(nodes_name)))
    z = np.zeros((len(nodes_name),len(nodes_name)))
    m = np.zeros((len(nodes_name),len(nodes_name)))
    for i in range(len(nodes_name)):
        for j in range(i+1, len(nodes_name)):
            if Dis[i,j] == 1 and i != j:
                if(any(gene_value[i,:])!=0 and any(gene_value[j,:])!=0):
                    P[i,j] = np.corrcoef(gene_value[i,:],gene_value[j,:])[0,1]
                else:
                    P[i,j] = 0
                temp = 0
                for k in range(len(nodes_name)):
                    if Dis[i,k] == 1 and Dis[j,k] ==1 and i!=j!=k:
                        temp+=1
                z[i,j] = temp
                m[i,j] = min(degree[i]-1,degree[j]-1)
                Ecc[i,j] = z[i,j]/m[i,j]
                Ecc[i,j] = Ecc[j,i] = 0.5*Ecc[i,j]+0.5*P[i,j]

    Ecc[np.isnan(Ecc)==1]=0
    lc = np.sum(Ecc,axis = 1)
    return lc 








