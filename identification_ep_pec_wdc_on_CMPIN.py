import pandas as pd
import numpy as np
import networkx as nx
import operator
from math import *
from operator import itemgetter, attrgetter
from math import *
from sklearn import preprocessing
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import auc
from sklearn.metrics import roc_auc_score
from PEC_method import pec
from WDC_method import wdc
from JDC_method import jdc
from UC_method import uc
species = 'yeast'
dataset = 'ydip'
pin = 'rdpin'
method = 'wdc'

Dis=np.load('model_detection/%s/%s/%s/CMPIN/final_dis.npy'%(species,dataset,pin))
nodes_name=np.load('model_detection/%s/%s/%s/CMPIN/final_nodes.npy'%(species,dataset,pin))
essential_nodes_name=np.load(r"model_detection/%s/%s/ess_gene_name.npy"%(species,dataset))
gene_value = np.load(r'model_detection/%s/%s/gene_value_profile.npy'%(species,dataset))

nodes_original=np.load("model_detection/%s/%s/nodes_name.npy"%(species,dataset))

f_dis = np.zeros((len(nodes_original),len(nodes_original)))
for i in range(len(Dis)):
    for j in range(i+1,len(Dis)):
        if Dis[i,j]==1:
            ind_1 = np.where(nodes_original == nodes_name[i])
            ind_2 = np.where(nodes_original == nodes_name[j])
            f_dis[ind_1[0], ind_2[0]] = f_dis[ind_2[0], ind_1[0]] = 1
            
lc = wdc(f_dis, nodes_original, gene_value)

y=np.zeros(len(nodes_original))
for i in range(len(nodes_original)):
    if nodes_original[i] in essential_nodes_name:
        y[i]=1

y_predict = np.array(lc)

#PRAUC
lr_precision, lr_recall, thresholds = precision_recall_curve(y, y_predict)
pr_auc = auc(lr_recall, lr_precision)
print('PR-AUC = %.4f' % (pr_auc))



ec_val=pd.DataFrame(lc,index=nodes_original)
aaa=ec_val.sort_values(by = 0,axis = 0,ascending = False)
edgelist_new=[]
edgelist_new=np.array(aaa.index)

edgelist_new600=[]
edgelist_new600=edgelist_new[:len(essential_nodes_name)]
essential = np.array(edgelist_new600)

all=[]
all100=[]
all200=[]
all300=[]
all400=[]
all500=[]
all600=[]
ess100=essential[:100]
ess200=essential[:200]
ess300=essential[:300]
ess400=essential[:400]
ess500=essential[:500]
ess600=essential[:600]
for id in essential_nodes_name:
    if id in ess100:
        all100.append(id)
print('top100:',len(all100))
for id in essential_nodes_name:
    if id in ess200:
        all200.append(id)
print('top200:',len(all200))
for id in essential_nodes_name:
    if id in ess300:
        all300.append(id)
print('top300:',len(all300))
for id in essential_nodes_name:
    if id in ess400:
        all400.append(id)
print('top400:',len(all400))
for id in essential_nodes_name:
    if id in ess500:
        all500.append(id)
print('top500:',len(all500))
for id in essential_nodes_name:
    if id in ess600:
        all600.append(id)
print('top600:',len(all600))
for id in essential:
    if id in essential_nodes_name:
        all.append(id)
print('top ess_num:',len(all))

