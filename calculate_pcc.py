import numpy as np
import pandas as pd

def pearson_model_and_orthologous(new_nodes, model_dete, orthologous_data):
    
    orthologous_name = np.array(orthologous_data.values[:,0])
    orthologous_val = np.array(orthologous_data.values[:,-1])
    orthologous_m = np.zeros(len(new_nodes))
    for i in range(len(new_nodes)):
        if new_nodes[i] in orthologous_name:
            ind = np.where(orthologous_name == new_nodes[i])
            orthologous_m[i] = orthologous_val[ind]
    
    model_m = np.zeros((len(new_nodes),len(model_dete)))# 模块的01矩阵
    
    for i in range(len(model_dete)):
        temp = model_dete[i]
        for j in range(len(temp)):
            ind_a=np.where(new_nodes == temp[j])
            model_m[ind_a,i]=1
        
        
    col = []
    col.extend(range(0,len(model_dete),1))
    col.append('ess')
    pear_frame = pd.DataFrame(data=np.hstack((model_m, orthologous_m.reshape(len(new_nodes),1))), columns = col)
    fianl_score = np.array(pear_frame.corr('pearson'))[-1,:-1]
    return fianl_score