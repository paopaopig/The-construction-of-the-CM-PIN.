import numpy as np
import networkx as nx


def max_connected_component(Dis, nodes_name, file):
    num_of_nodes = len(nodes_name)
    G = nx.Graph()
    G.add_nodes_from(nodes_name) 
    for i in range(num_of_nodes):
        for j in range(num_of_nodes):
            if(Dis[i][j]==1):
                G.add_edge(nodes_name[i],nodes_name[j])
                
    #print("The number of connected components is：",nx.number_connected_components(G))
    C = sorted(nx.connected_components(G), key=len, reverse=True)
    a=[]
    final_nodes = []
    for i in range(len(nodes_name)):#原始节点的数量
        if nodes_name[i] not in list(C[0]):
            a.append(i)
        if nodes_name[i] in list(C[0]):
            final_nodes.append(nodes_name[i])
                
    Dis = np.delete(Dis, a, axis=1)
    final_Dis = np.delete(Dis, a, axis=0) 
    
    print("The number of interactions in the maximum connected subgraph：",sum(sum(final_Dis))/2)
    print("The number of nodes in the maximum connected subgraph：：",len(final_nodes))
    
    f = open(file,'a')
    for i in range(len(final_nodes)):
        for j in range(i+1,len(final_nodes)):
            if final_Dis[i,j]==1:
                f.writelines(final_nodes[i])
                f.writelines('\t')
                f.writelines(final_nodes[j])
                f.writelines('\t')
                f.writelines(str(1))
                f.writelines('\n')
    f.close()
    print("The extraction of the maximum connected subgraph is complete!")
    return final_Dis, final_nodes

#if __name__ == "__main__":
    #q, p, path = max_connected_component(Dis, nodes_name, 'ydip_spin')