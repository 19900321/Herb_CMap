import pandas as pd
from cal_proximity_shortest import closest_shortest_matrix

def cal_LR_disase_gene_DEG(G, deg_nodes, disease_nodes):

    # calculate taregt-DEG distance
    DEG_dis_gene_distance_matrix = closest_shortest_matrix(G, deg_nodes, disease_nodes, weight=True)
    
    # Apply Local radility alogrithm  base on clsoet
    DEG_dis_gene_lr_score = pd.DataFrame.from_dict(DEG_dis_gene_distance_matrix.sum(axis=1).to_dict(),orient='index', columns=['LR'] )
    DEG_dis_gene_lr_score.insert(0,'Disease-related_gene',  DEG_dis_gene_lr_score.index)
    
    return DEG_dis_gene_distance_matrix, DEG_dis_gene_lr_score
