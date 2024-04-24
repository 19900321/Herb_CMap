import networkx as nx
import pandas as pd
import tqdm

# calcualte shortest distance
def closest_shortest_matrix(G, A, B, weight=True):
    dis_pd = pd.DataFrame(columns=A, index=B)
    for node1 in tqdm.tqdm(A):
        for node2 in B:
            print('1')
            try:
                if weight==True:
                    # distance, weighted_shortest_path = cauculate_four_step_closet(G, node1, node2)
                    #
                    # '''
                    # the following methods is too time-consuming
                    # '''
                    # #
                    distance = nx.shortest_path_length(G, node1, node2, weight='w')
                else:
                    distance = nx.shortest_path_length(G, node1, node2)
            except nx.NetworkXNoPath:
                distance = None
            dis_pd.loc[node2, node1] = distance
    return dis_pd
