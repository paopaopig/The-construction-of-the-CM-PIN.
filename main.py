import numpy as np
import pandas as pd
from refined_by_max_con import max_connected_component
from fast_unfolding import loadData, fast_unfolding
from calculate_pcc import pearson_model_and_orthologous
from calculate_nsl import model_sub_score
from calculate_topo import model_tp_score
import os

species = 'human'
dataset = 'hdip'
pin = 'rdpin'

Dis = np.load(r"model_detection/%s/%s/%s/%s.npy"%(species,dataset,pin,pin))
nodes_name = np.load(r"model_detection/%s/%s/nodes_name.npy"%(species,dataset))
essential_nodes_name=np.load("model_detection/%s/%s/ess_gene_name.npy"%(species,dataset))

#********The first step: get the maximally connected subgraph and its new node set and adjacency matrix********
file_txt = 'model_detection/%s/%s/%s/%s_%s.txt'%(species,dataset,pin,dataset, pin)
# Write the maximum connected subgraph to a txt file
if os.access(file_txt, os.F_OK):
    new_dis = np.load('model_detection/%s/%s/%s/new_dis_%s_%s.npy'%(species,dataset,pin,dataset,pin), allow_pickle=True)
    new_nodes = np.load('model_detection/%s/%s/%s/new_nodes_%s_%s.npy'%(species,dataset,pin,dataset,pin), allow_pickle=True)
else:
    new_dis, new_nodes = max_connected_component(Dis, nodes_name, file_txt)
    np.save('model_detection/%s/%s/%s/new_dis_%s_%s.npy'%(species,dataset,pin,dataset,pin),new_dis)
    np.save('model_detection/%s/%s/%s/new_nodes_%s_%s.npy'%(species,dataset,pin,dataset,pin),new_nodes)
print("Complete the first step to extract the maximum connected subgraph！")

#*******The second step: the maximal connected subgraph is divided into modules********
file_model_dete = 'model_detection/%s/%s/%s/model_detection_%s_%s.npy'%(species,dataset,pin,dataset,pin)
if os.access(file_model_dete, os.F_OK):
    model_dete = np.load('model_detection/%s/%s/%s/model_detection_%s_%s.npy'%(species,dataset,pin,dataset,pin), allow_pickle=True)
else:
    vector_dict, edge_dict = loadData(file_txt)
    model_dete = fast_unfolding(vector_dict, edge_dict)
    np.save('model_detection/%s/%s/%s/model_detection_%s_%s.npy'%(species,dataset,pin,dataset,pin),model_dete)
print("Complete the second step, the division of modules！")
#********Step 3: Score each module (three scores)********
# 1.Pearson correlation coefficients were obtained for each module and protein homology score
orthologous_data = pd.read_csv('model_detection/%s/%s/orthologous.csv'%(species,dataset),header=None)
pcc = pearson_model_and_orthologous(new_nodes, model_dete, orthologous_data).reshape(len(model_dete),1)
print("The pcc is complete！")
# 2.The subcellular localization score of each module was obtained
sub_data_path = 'model_detection/%s/subcellular.txt'%species
nsl = model_sub_score(new_nodes, model_dete, sub_data_path).reshape(len(model_dete),1)
print("The nsl is complete！")
# 3.The topological feature score of each module is obtained
topo = model_tp_score(new_dis, new_nodes, model_dete).reshape(len(model_dete),1)
print("The topo is complete！")
# Step 4: Select the modules you want to keep
ind = np.arange(len(model_dete))
final_df = pd.DataFrame(data = np.hstack((pcc,nsl,topo)), index = ind, columns = ['pcc','nsl','topo'])

# Set three thresholds
th1 = 0.02
th2 = 0.4
th3 = 0
df_pcc = final_df[(final_df["pcc"] >= th1)]
df_nsl = final_df[(final_df["nsl"] >= th2)]
df_topo = final_df[(final_df["topo"] <= th3)]

final_selected_model = set(df_nsl.index).difference(set(df_topo.index)).union(set(df_pcc.index))
final_deleted_model = list(set(ind).difference(final_selected_model))


# Step 5: Delete nodes in final_deteted_model
del_nodes = []
for i in range(len(final_deleted_model)):
    temp = model_dete[final_deleted_model[i]]
    del_nodes.extend(temp)

del_nodes_ind = []
final_nodes = []
for i in range(len(new_nodes)):
    if new_nodes[i] in del_nodes:
        del_nodes_ind.append(i)
    if new_nodes[i] not in del_nodes:
        final_nodes.append(new_nodes[i])

new_dis= np.delete(new_dis, del_nodes_ind, axis=1)
final_Dis=np.delete(new_dis, del_nodes_ind, axis=0)

np.save('model_detection/%s/%s/%s/CMPIN/final_nodes.npy'%(species,dataset,pin),final_nodes)
np.save('model_detection/%s/%s/%s/CMPIN/final_dis.npy'%(species,dataset,pin),final_Dis)




    