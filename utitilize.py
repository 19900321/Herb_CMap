
import pandas as pd

def pre_suhuang_herb_target_expression():
    ingre_target = pd.read_csv('source_data/weighted_rat_ppi/SH_ingre_tar_withscore.csv')

    # map target resource from ingre to herb to property
    # read herb property
    herb_anno = pd.read_excel('source_data/weighted_rat_ppi/SUHUANG_h_property.xlsx', sheet_name='Sheet1')
    
    # read herb_ingre
    herb_ingre = pd.read_csv('source_data/weighted_rat_ppi/SUHUANG_h_ingre.csv')
    
    # map pro to ingredient
    h_i_anno = pd.merge(herb_ingre, herb_anno, right_on='Chinese Name', left_on='中文名', how='left')

    # map property to ingre target
    ingre_h_i_anno = pd.merge(ingre_target, h_i_anno, left_on='inchikey', right_on='inchikey', how='inner')
    ingre_h_i_anno = ingre_h_i_anno.drop_duplicates()

    # map gene expression
    gene_expre_pd_all = pd.read_csv('source_data/weighted_rat_ppi/SHVSM_expr_result_all.csv')
    ingre_h_i_anno_expression = pd.merge(ingre_h_i_anno, gene_expre_pd_all, left_on='ENSEMBL', right_on='ENSEMBL',
                                         how='inner')
    ingre_h_i_anno_expression.to_csv('source_data/weighted_rat_ppi/SUHUANG_h_ingre_t_expression.csv', index=None)


def prepare_target_anno(target_threshold, ingre_anno_path, nodes_ppi, save_path):
    # read ingre_taregt dile
    ingre_h_i_anno_expression = pd.read_csv(ingre_anno_path)

    # select how many ingredient to select
    ingre_target_selected = ingre_h_i_anno_expression.loc[ingre_h_i_anno_expression['combined_score'] >= target_threshold,]
    ingre_with_target = list(ingre_target_selected['inchikey'].unique())
    print('we get {} target in total'.format(len(ingre_with_target)))# 245
    target_nodes = list(ingre_target_selected['ENSEMBL'].unique())
    target_nodes_2 = [t for t in target_nodes if t in nodes_ppi]


    # Use pivot_table() to create the matrix
    target_herb_matrix = ingre_target_selected[['ENSEMBL', 'Chinese Name']].pivot_table(index='ENSEMBL', columns= 'Chinese Name', aggfunc=len, fill_value=0)

    target_siqi_matrix = ingre_target_selected[['ENSEMBL', 'Siqi']].pivot_table(index='ENSEMBL', columns='Siqi', aggfunc=len, fill_value=0)
    target_wuwei_matrix = ingre_target_selected[['ENSEMBL', 'Wuwei']].pivot_table(index='ENSEMBL', columns='Wuwei', aggfunc=len,
                                                                      fill_value=0)
    target_meridians_matrix = ingre_target_selected[['Meridians', 'ENSEMBL']].pivot_table(index='ENSEMBL', columns='Meridians', aggfunc=len,
                                                                        fill_value=0)

    # SAVE OUT
    ingre_target_selected.to_csv('{}/target_info.csv'.format(save_path), index=None)
    target_herb_matrix.to_csv('{}/target_herb_matrix.csv'.format(save_path))
    target_siqi_matrix.to_csv('{}/target_siqi_matrix.csv'.format(save_path))
    target_wuwei_matrix.to_csv('{}/target_wuwei_matrix.csv'.format(save_path))
    target_meridians_matrix.to_csv('{}/target_meridians_matrix.csv'.format(save_path))
    return target_nodes_2


def prepare_all_DEGs(deg_path, nodes_ppi):
    # load ingre_target， DEGs
    deg_pd = pd.read_csv(deg_path)
    deg_nodes = list(deg_pd['ENSEMBL'].unique())
    deg_nodes_2 = [t for t in deg_nodes if t in nodes_ppi]
    return deg_nodes_2


def prepare_cva_disease_gene_(disease_gene_path, nodes_ppi):
    disease_related_gene = pd.read_csv(disease_gene_path)
    dissease_nodes = list(disease_related_gene['ENSEMBL'].unique())
    disease_nodes_2 = [t for t in dissease_nodes if t in nodes_ppi]
    return disease_nodes_2


# only keep DEGs in disease and suhuang overalp 332, and in five pathways
def prepare_pathway_genes_anno(nodes_ppi, pathway_selelcted_path, gene_expression_path):
    pathways = pd.read_excel(pathway_selelcted_path, sheet_name='Sheet1')
    pathways['geneID_list_list'] = pathways['geneID'].apply(lambda x:x.split('/'))

    # Use get_dummies() to create new columns with 1 and 0
    pathways_gene_matrix = pd.get_dummies(pathways['geneID_list_list'].apply(pd.Series).stack()).sum(level=0)

    # Concatenate the original dataframe and the new columns
    pathways_gene_matrix.index = pathways['Description']
    pathways_gene_matrix = pathways_gene_matrix.T
    pathways_gene_matrix['SYMBOL'] = pathways_gene_matrix.index

    # read the gene expression value
    gene_expre_pd = pd.read_csv(gene_expression_path)
    gene_expre_pd = gene_expre_pd.sort_values(by='log2FoldChange', ascending=True)
    gene_expre_pd.to_csv('result/suhuang/target_gene_heat/gene_info_all.csv', index= None)

    # rename the gene to ens in pathways_expand
    # GENERATE GENE_SYMBL_ENS DICT
    gene_info = pd.merge(gene_expre_pd, pathways_gene_matrix, right_on='SYMBOL', left_on='SYMBOL', how='inner')
    gene_info = gene_info[gene_info['ENSEMBL'].isin(nodes_ppi)]
    gene_info = gene_info.sort_values(by='log2FoldChange', ascending=True)

    # replace symbol to ensembl
    pathways_gene_matrix = gene_info[['ENSEMBL', 'baseMean', 'log2FoldChange', 'lfcSE',
       'stat', 'pvalue', 'padj', 'group', 'Rheumatoid arthritis',
       'IL-17 signaling pathway', 'TNF signaling pathway',
       'Toll-like receptor signaling pathway', 'Asthma',
       'T cell receptor signaling pathway',
       'PD-L1 expression and PD-1 checkpoint pathway in cancer']]
    return gene_info, pathways_gene_matrix




  