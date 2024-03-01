import numpy as np
import pandas as pd

def model_tp_score(new_dis, new_nodes, model_dete):
    model_num = len(model_dete)
    new_model_0_1 = []
    for i in range(model_num):
        #print(i)
        temp = model_dete[i].copy()
        for j in range(len(temp)):
            for k in range(len(new_nodes)):
                if temp[j]==new_nodes[k]:
                    temp[j]=k
                
        new_model_0_1.append(temp)        
    
    Dis_new=np.zeros((model_num,model_num))

    for i in range(model_num):
        for j in range(i+1,model_num):
            a = new_model_0_1[i]
            b = new_model_0_1[j]
            for m in range(len(a)):
                for n in range(len(b)):
                    Dis_new[i,j]+=new_dis[a[m],b[n]]
    for i in range(model_num):
        for j in range(i+1,model_num):
            Dis_new[j,i]=Dis_new[i,j]
            

    community_nodes=np.zeros(model_num)
    for i in range(model_num):
        temp=new_model_0_1[i]
        community_nodes[i]=len(temp)
    community_nodes=community_nodes
    
    Dis_new2=np.zeros(model_num)

    for i in range(model_num):
        a=new_model_0_1[i]        
        for m in range(len(a)):
            for n in range(m+1,len(a)):
                Dis_new2[i]+=new_dis[a[m],a[n]]
                
    final=np.zeros(model_num)
    for i in range(model_num):
        t=0
        for j in range(model_num):
            t+=Dis_new[i,j]
        final[i]=(Dis_new2[i]-t)/community_nodes[i]
            
    return final

                