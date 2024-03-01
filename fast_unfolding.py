import numpy as np

def loadData(filePath):
    f = open(filePath)
    vector_dict = {}
    edge_dict = {}
    for line in f.readlines():
        lines = line.strip().split("\t")
        for i in range(2):
            if lines[i] not in vector_dict:
                vector_dict[lines[i]] = True
                edge_list = []
            else:
                edge_list = edge_dict[lines[i]]
            edge_list.append(lines[1-i]+":"+lines[2])
            edge_dict[lines[i]] = edge_list
    return vector_dict, edge_dict




def modularity(vector_dict, edge_dict):
    Q = 0.0
    m = 0
    for i in edge_dict.keys():
        edge_list = edge_dict[i]
        for j in range(len(edge_list)):
            l = edge_list[j].strip().split(":")
            m += float(l[1].strip())
    community_dict = {}
    for i in vector_dict.keys():
        if vector_dict[i] not in community_dict:
            community_list = []
        else:
            community_list = community_dict[vector_dict[i]]
        community_list.append(i)
        community_dict[vector_dict[i]] = community_list
    for i in community_dict.keys():
        sum_in = 0.0
        sum_tot = 0.0
        vector_list = community_dict[i]
        for j in range(len(vector_list)):
            link_list = edge_dict[vector_list[j]]
            tmp_dict = {}
            for link_mem in link_list:
                l = link_mem.strip().split(":")
                tmp_dict[l[0]] = l[1]
            for k in range(0, len(vector_list)):
                if vector_list[k] in tmp_dict:
                    sum_in += float(tmp_dict[vector_list[k]])
        for vec in vector_list:
            link_list = edge_dict[vec]
            for i in link_list:
                l = i.strip().split(":")
                sum_tot += float(l[1])
        Q += ((sum_in / m) - (sum_tot/m)*(sum_tot/m))
    return Q



def change_community(vector_dict, edge_dict, Q):
    vector_tmp_dict = {}
    for key in vector_dict:
        vector_tmp_dict[key] = vector_dict[key]
    for key in vector_tmp_dict.keys():
        neighbor_vector_list = edge_dict[key]
        for vec in neighbor_vector_list:
            ori_com = vector_tmp_dict[key]
            vec_v = vec.strip().split(":")
            if ori_com != vector_tmp_dict[vec_v[0]]:
                vector_tmp_dict[key] = vector_tmp_dict[vec_v[0]]
                Q_new = modularity(vector_tmp_dict, edge_dict)
                if (Q_new - Q) > 0:
                    Q = Q_new
                else:
                    vector_tmp_dict[key] = ori_com
    return vector_tmp_dict, Q



def modify_community(vector_dict):
    community_dict = {}
    community_num = 0
    for community_values in vector_dict.values():
        if community_values not in community_dict:
            community_dict[community_values] = community_num
            community_num += 1
    for key in vector_dict.keys():
        vector_dict[key] = community_dict[vector_dict[key]]
    return community_num



def rebuild_graph(vector_dict, edge_dict, community_num):
    vector_new_dict = {}
    edge_new_dict = {}
    community_dict = {}
    for key in vector_dict.keys():
        if vector_dict[key] not in community_dict:
            community_list = []
        else:
            community_list = community_dict[vector_dict[key]]
        community_list.append(key)
        community_dict[vector_dict[key]] = community_list
    for key in community_dict.keys():
        vector_new_dict[str(key)] = str(key)
    for i in community_dict.keys():
        sum_in = 0.0
        vector_list = community_dict[i]
        if '0' in vector_list:
            print(vector_list)
            print(i)
        for j in range(0,len(vector_list)):
            link_list = edge_dict[vector_list[j]]
            tmp_dict = {}
            for link_mem in link_list:
                l = link_mem.strip().split(":")
                tmp_dict[l[0]] = l[1]
            for k in range(0, len(vector_list)):
                if vector_list[k] in tmp_dict:
                    sum_in += float(tmp_dict[vector_list[k]])
        inner_list = []
        inner_list.append(str(i) + ":" + str(sum_in))
        edge_new_dict[str(i)] = inner_list
    community_list = list(community_dict.keys())
    for i in range(len(community_list)):
        for j in range(len(community_list)):
            if i != j:
                sum_outer = 0.0
                member_list_1 = community_dict[community_list[i]]
                member_list_2 = community_dict[community_list[j]]
                for i_1 in range(len(member_list_1)):
                    tmp_dict = {}
                    tmp_list = edge_dict[member_list_1[i_1]]
                    for k in range(len(tmp_list)):
                        tmp = tmp_list[k].strip().split(":")
                        tmp_dict[tmp[0]] = tmp[1]
                    for j_1 in range(len(member_list_2)):
                        if member_list_2[j_1] in tmp_dict:
                            sum_outer += float(tmp_dict[member_list_2[j_1]])
                if sum_outer != 0:
                    inner_list = edge_new_dict[str(community_list[i])]
                    inner_list.append(str(j) + ":" + str(sum_outer))
                    edge_new_dict[str(community_list[i])] = inner_list
    return vector_new_dict, edge_new_dict, community_dict



def fast_unfolding(vector_dict, edge_dict):
    for i in vector_dict.keys():
        vector_dict[i] = i
    Q = modularity(vector_dict, edge_dict)
    Q_new = 0.0
    while (Q_new != Q):
        Q_new = Q
        vector_dict, Q = change_community(vector_dict, edge_dict, Q)
    community_num = modify_community(vector_dict)
    print("Q = ", Q)
    print(community_num)
    '''
    print("vector_dict.key : ", vector_dict.keys())
    print("vector_dict.value : ", vector_dict.values())
    '''
    Q_best = Q
    while True:
        print("\n rebuild")
        vector_new_dict, edge_new_dict, community_dict = rebuild_graph(vector_dict, edge_dict, community_num)
        print("community_dict : ", community_dict)
        
        Q_new = 0.0
        while (Q_new != Q):
            Q_new = Q
            vector_new_dict, Q = change_community(vector_new_dict, edge_new_dict, Q)
        community_num = modify_community(vector_new_dict)

        print("Q = ", Q)
        print("community_num : ", community_num)
        if (Q_best == Q):
            break
        Q_best = Q

        vector_result = {}
        for key in community_dict.keys():
            value_of_vector = community_dict[key]
            for i in range(len(value_of_vector)):

                vector_result[value_of_vector[i]] = str(vector_new_dict[str(key)])
        for key in vector_result.keys():
            vector_dict[key] = vector_result[key]

    vector_result = {}
    for key in community_dict.keys():
        value_of_vector = community_dict[key]
        for i in range(len(value_of_vector)):
            vector_result[value_of_vector[i]] = str(vector_new_dict[str(key)])
    for key in vector_result.keys():
        vector_dict[key] = vector_result[key]
    print("Q_best : ", Q_best)

    community_dict = np.array(community_dict).tolist()
    community = list(community_dict.values())
    print("模块划分完成！")
    return community

