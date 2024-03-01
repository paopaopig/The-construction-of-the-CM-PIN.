import numpy as np
import pandas as pd

def model_sub_score(new_nodes,model_dete,sub_data_path):
    sub = pd.read_csv(sub_data_path, sep='\t',header=None)
    sub_num = sum(1 for line in open(sub_data_path))
    print(sub_num)
    sub_name=np.array(sub.iloc[:,0]).reshape(sub_num,1).ravel()
    subcell_11=np.array(sub.iloc[:,1]).reshape(sub_num,1).ravel()

    
    sub_nodes=[]
    sub_room=[]
    for i in range(sub_num):
        if sub_name[i] in new_nodes:
            sub_nodes.append(sub_name[i])
            sub_room.append(subcell_11[i])
    
    sub_nodes = np.array(sub_nodes)
    sub_room = np.array(sub_room)
    
    final_score = np.zeros(len(model_dete))
    for i in range(len(model_dete)):
        temp_comm = model_dete[i]     
        for j in range(len(temp_comm)):
            ind_p = np.array(np.where(sub_nodes==temp_comm[j])).ravel()
            for k in range(len(ind_p)):
                if sub_room[ind_p[k]]=='Nucleus':
                    final_score[i]+=1
        final_score[i] = final_score[i]/len(temp_comm)
    
    return final_score