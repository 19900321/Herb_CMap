from scipy import stats

def pre_weighted_ppi(ppi_pd, suhuang_expressionp_pd):
    
    # prepare ensemble weight dictionary
    weight_dict = dict(zip(suhuang_expressionp_pd['ENSEMBL'], abs(suhuang_expressionp_pd['log2FoldChange'])))
    
    # apply weight to source and target to get averaged value
    ppi_pd['source_gene_expre_value'] = ppi_pd['source_gene'].apply(lambda x:weight_dict.get(x))
    ppi_pd['target_gene_expre_value'] = ppi_pd['target_gene'].apply(lambda x:weight_dict.get(x))
    ppi_pd['geneexpre_value'] = (ppi_pd['source_gene_expre_value'] + ppi_pd['target_gene_expre_value'])/2

    # transform data from non-normally distributed dataset into a more normally distributed one.
    '''
    A box-cox transformation is a commonly used method for transforming a non-normally distributed dataset into a more normally distributed one.
    '''
    fitted_data, _ = stats.boxcox(ppi_pd['geneexpre_value'])

    # use max -value to return the expression to distance weight
    ppi_pd['geneexpre_value'] = -fitted_data+fitted_data.max()
    
    # select the sorce weight and target as columns
    ppi_pd_selected = ppi_pd[['source_gene','geneexpre_value','target_gene']]
    
    return ppi_pd_selected
