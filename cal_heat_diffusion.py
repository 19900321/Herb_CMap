import pandas as pd
from cal_proximity_shortest import closest_shortest_matrix
import numpy as np
from scipy.linalg import expm
import os

def cal__DEG_target_distance_matrix(G, deg_nodes, target_nodes, taregt_DEG_path):

    # calculate taregt-DEG distance
    if os.path.exists(taregt_DEG_path):
        DEG_target_distance_matrix = pd.read_csv(taregt_DEG_path, index_col = 0)
    else:
        DEG_target_distance_matrix = closest_shortest_matrix(G, deg_nodes, target_nodes, weight=True)
        DEG_target_distance_matrix.to_csv(taregt_DEG_path)
    
    # Rank all values
    sorted_values  = sorted([i for i in DEG_target_distance_matrix.to_numpy().flatten() if isinstance(i, float)]) 
    
    # found the top 5% value for define threshold
    threshold = sorted_values[int(len(sorted_values)*0.05)]
    print("cut threshold is {}".format(threshold))
    
    # Apply the condition,  map 5% target-gene as 1 and the others as 0
    DEG_target_distance_matrix_binary = (DEG_target_distance_matrix <= threshold).astype(int)
    
    # mask the weight 
    DEG_target_distance_matrix_binary = DEG_target_distance_matrix_binary*DEG_target_distance_matrix

    # cut unrelated nodes
    DEG_target_distance_matrix_binary = DEG_target_distance_matrix_binary.loc[DEG_target_distance_matrix_binary.sum(axis=1)!=0, DEG_target_distance_matrix_binary.sum(axis=0)!=0]
    
    return DEG_target_distance_matrix, DEG_target_distance_matrix_binary
    
    
def pre_pare_DEG_plus_target_matrix(DEG_target_distance_matrix_binary):
    DEG_Target_names = list(DEG_target_distance_matrix_binary.columns) + list(DEG_target_distance_matrix_binary.index)
    DEG_target_all_matrix = pd.DataFrame(np.diag([1]*len(DEG_Target_names)))
    DEG_target_all_matrix.columns = DEG_Target_names
    DEG_target_all_matrix.index =  DEG_Target_names
    
    for gene1 in DEG_target_distance_matrix_binary.index:
        for gene2 in DEG_target_distance_matrix_binary.columns:
            DEG_target_all_matrix.loc[gene1, gene2] = DEG_target_distance_matrix_binary.loc[gene1, gene2]
            DEG_target_all_matrix.loc[gene2, gene1] = DEG_target_distance_matrix_binary.loc[gene1, gene2]
    
    return DEG_target_all_matrix


def heat_diffusion(DEG_target_all_matrix, gene_expression_h, t):        
    # Apply heat diffusion alogrithm  base on shortest distance

    # Compute the diagonal matrix E where each element is the number of connections (degree) of a node in C
    degree = np.sum(DEG_target_all_matrix, axis=0)
    E = np.diag(degree)
    L = E - DEG_target_all_matrix

    # Diffusion process with convergence 
    epsilon = 10  # Convergence threshold
    max_iter = 3  # Maximum number of iterations to prevent infinite loops
    
    # Initial heat vector h representing the gene expression change of targets and DEGs
    current_d = gene_expression_h.copy()

    for i in range(max_iter):
        next_d = current_d @ expm(-L * t)
        #print('The loss is ', np.linalg.norm(next_d - current_d))
        # Check if the system has converged
        if np.linalg.norm(next_d - current_d) < epsilon:
            break
        current_d = next_d
        #print('cuurent d', current_d)

    print("Converged after {} iterations.".format(i+1))
    print("The final diffusion vector is:", current_d)
    
    # prepare as dataframe
    diffusion_result_pd = pd.DataFrame({'genes':DEG_target_all_matrix.index,
                                        'diffusion_input': gene_expression_h,
                                        'diffusion_output_heat':current_d})
    return  diffusion_result_pd


def cal_target_DEG_heat_score(G, suhuang_expressionp_pd, deg_nodes, target_nodes, t, taregt_DEG_path):
    ''' # DEG_target_distance_pd = pd.read_csv('result/suhuang/target_gene_dis/target_DEG_distance_z_p.csv')
    # DEG_target_distance_matrix = DEG_target_distance_pd.pivot_table(columns='DEG',index='target', values= 'distance')
    # taregt_DEG_path  = 'result/suhuang/target_gene_dis/target_gene_dis_all.csv'
    # DEG_target_distance_matrix.to_csv(taregt_DEG_path)'''
   
    # calculate distance
    DEG_tar_distance_matrix, DEG_tar_dist_mat_binary = cal__DEG_target_distance_matrix(G, deg_nodes, target_nodes, taregt_DEG_path)
    
    # prepare the large connected target-DEG network with distance as 
    DEG_tar_all_matrix = pre_pare_DEG_plus_target_matrix(DEG_tar_dist_mat_binary)
    
    # select the gene expression 
    gene_expression_h = np.array(abs(suhuang_expressionp_pd.loc[DEG_tar_all_matrix.columns,'log2FoldChange']))
    
    diffusion_result_pd = heat_diffusion(DEG_tar_all_matrix, gene_expression_h, t)
  
    return DEG_tar_distance_matrix, DEG_tar_dist_mat_binary, diffusion_result_pd