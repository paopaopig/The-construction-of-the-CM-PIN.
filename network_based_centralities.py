import pandas as pd
import numpy as np
import networkx as nx
from math import *

def lid(Dis,nodes):    
    lids={}
    lid=0
    for i in range(len(nodes)):
        b=np.argwhere(Dis[i,:]==1)[:,0]
        if len(b):
            C=Dis[np.ix_(b,b)]
            Cv=np.sum(C)/2
            C_count=np.sum(C,axis=0)
            En=(list(C_count).count(0))
            Enb=len(C_count)-En
            if Enb==0:
                lid=0
            else:
                lid=Cv/Enb           
            lids.update({nodes[i]: lid})
        else:
            lids.update({nodes[i]: 0})
    delta_df = pd.DataFrame(columns=['顶点蛋白质id', 'LID值'], index=None)
    for key, v in enumerate(lids):
        delta_df.loc[key] = [v, lids[v]]
    return np.array(delta_df['LID值'])

def lac(Dis,nodes):
    lacvalue=[0 for i in range(len(nodes))]
    for i in range(len(nodes)):
        linju=np.argwhere(Dis[i,:]==1)[:,0]
        if len(linju):
            C=Dis[np.ix_(linju,linju)]
            Cv=np.sum(C)
            lacvalue[i]=Cv/sum(Dis[i,:])
    return lacvalue

def dmnc(Dis,nodes):
    dmncvalue=[0 for i in range(len(nodes))]
    for i in range(len(nodes)):   
        linju=np.argwhere(Dis[i,:]==1)[:,0]
        if len(linju):
            C=Dis[np.ix_(linju,linju)]
            Eci=np.sum(C)/2
            Nci=sum(Dis[i,:])**1.7
            dmncvalue[i]=Eci/Nci
    return dmncvalue


def nc(Dis,nodes):
    uu=[0 for i in range(1)]
    ncvalue=[0 for i in range(len(nodes))]
    for i in range(len(nodes)):
        linju=np.argwhere(Dis[i,:]==1)[:,0]
        nccount=0
        Zuv=0
        if len(linju)>1:
            for j in range(len(linju)):
                uu.clear
                Ecc=0
                u=linju[j]
                uu[0]=u
                C=Dis[np.ix_(uu,linju)]
                Zuv=np.sum(C)
                if Zuv>0:
                    Ecc=Zuv/min(sum(Dis[u,:])-1,sum(Dis[i,:])-1)
                nccount=nccount+Ecc
            ncvalue[i]=nccount
        else:
            ncvalue[i]=0
    return ncvalue

def tp(Dis,nodes):
    G = nx.Graph()
    G.add_nodes_from(nodes)  
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if(Dis[i][j]==1):
                G.add_edge(nodes[i],nodes[j])
    tpvalue=[0 for i in range(len(nodes))]
    for i in range(len(nodes)):
        TPcount=0
        if(i%2000==0):
            print(i)
        for j in range(len(nodes)):
            TPj=0
            if nx.has_path(G,nodes[i],nodes[j]):
                P=nx.shortest_path_length(G,source=nodes[i],target=nodes[j])               
                P=P/0.9428
                P=P**2
                P=P*(-1)
                TPj=pow(e,P)
                TPcount=TPcount+TPj              
        tpvalue[i]=TPcount    
    return tpvalue


def dc(Dis,nodes):
    G = nx.Graph()
    G.add_nodes_from(nodes)  
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if(Dis[i][j]==1):
                G.add_edge(nodes[i],nodes[j])  
    dc=[0 for i in range(len(nodes))]
    dc = nx.algorithms.centrality.degree_centrality(G)
    dc_val=[]
    dc_val=dc.values()
    dc_val=list(dc_val)
    dc_val=np.array(dc_val)
    return dc_val


def cc(Dis,nodes):
    G = nx.Graph()
    G.add_nodes_from(nodes)  
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if(Dis[i][j]==1):
                G.add_edge(nodes[i],nodes[j])
    cc=[0 for i in range(len(nodes))]
    cc=nx.closeness_centrality(G)
    cc_val=[]
    cc_val=cc.values()
    cc_val=list(cc_val)
    cc_val=np.array(cc_val)
    return cc_val

def bc(Dis,nodes):
    G = nx.Graph()
    G.add_nodes_from(nodes) 
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if(Dis[i][j]==1):
                G.add_edge(nodes[i],nodes[j])
    bc=[0 for i in range(len(nodes))]  
    bc=nx.betweenness_centrality(G)
    bc_val=[]
    bc_val=bc.values()
    bc_val=list(bc_val)
    bc_val=np.array(bc_val)
    return bc_val

def pr(Dis,nodes):
    G = nx.Graph()
    G.add_nodes_from(nodes) 
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if(Dis[i][j]==1):
                G.add_edge(nodes[i],nodes[j])
    pr=[0 for i in range(4746)]
    pr=nx.pagerank(G, alpha=0.85)
    pr_val=[]
    pr_val=pr.values()
    pr_val=list(pr_val)
    pr_val=np.array(pr_val)
    return pr_val

def lr(Dis,nodes):
    activate_ppinet1=np.hstack((np.ones(shape=(len(nodes),1)),Dis))
    activate_ppinet1=np.vstack((np.ones(shape=(1,len(nodes)+1)),activate_ppinet1))
    tempLR=[0]*(len(nodes)+1)
    LR=[1]*(len(nodes)+1)
    LR[0]=0
    while True:
        for i in range(len(nodes)+1):
            s=0
            for j in range(i+1,len(nodes)+1,1):
                if activate_ppinet1[i,j]==1:
                    s+=((1/sum(activate_ppinet1[j,:]))*LR[j])
            tempLR[i]=s
        if tempLR==LR:
            
            break
        LR=tempLR
    avg=LR[0]/len(nodes)
    LR.pop(0)
    LR=np.array(LR)+avg
    result=LR.reshape(len(nodes),1).astype('float64')

    return result






    