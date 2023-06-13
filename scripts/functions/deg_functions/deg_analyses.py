import scipy
import warnings
import anndata2ri
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sb
import decoupler as dc
from scipy import sparse
from anndata import AnnData
from tabnanny import verbose
import matplotlib.pyplot as plt
from gsva_prep import prep_gsva
from functions import pathway_analyses
from typing import Optional, Union
from matplotlib.pyplot import rcParams
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import train_test_split
from rpy2.robjects.conversion import localconverter


def prep_anndata(adata: AnnData, map_meta: bool=True,
                metadata: str = '../data/raw/mathys_pfc/mathys_pfc_metadata.csv',
                reference_group: str = 'no', test_groups: list = ['late', "early"],
                subject_id: str = 'Subject',
                grouped_test: bool = True):
    """
    Add a disease group annotation to AnnData based on subject pathology status and groups.

    Parameters:
    -----------
    adata : AnnData
        AnnData object to be annotated.
    map_meta: bool, optional,
        Wether to map subjects/cells to pathology.group provided in metadata.
    metadata : str, optional
        Path to metadata file containing subject pathology group labels. Default is '../data/raw/mathys_pfc/mathys_pfc_metadata.csv'.
    reference_group : str, optional
        Name of the reference group for comparison. Default is 'no'.
    test_groups : list, optional
        List of pathology group names to be compared against reference group. Default is ['late', 'early'].
    subject_id : str, optional
        Name of the column containing subject IDs in the metadata file. Default is 'Subject'.
    grouped_test: bool, optional
        Whether to group tested levels into one group.

    Returns:
    --------
    adata : AnnData
        The annotated AnnData object.
    """

    # Check if pathology.group column is present in the adata object
    if map_meta==True:
        # Read metadata file and extract subject ID and pathology group labels
        meta = pd.read_csv(metadata)
        meta = meta.astype(str)
        meta = dict(zip(meta[subject_id], meta['pathology.group']))
        
        # Map the pathology group labels to the adata object using subject IDs
        adata.obs['pathology.group'] = adata.obs[subject_id].map(meta)

    # Add a disease group annotation to the adata object based on reference and test groups
    adata.obs['disease_group'] = None
    adata.obs.loc[adata.obs['pathology.group'] == reference_group, 'disease_group'] = reference_group+'-pathology'

    if grouped_test:
        adata.obs.loc[adata.obs['pathology.group'].isin(test_groups), 'disease_group'] = 'AD-pathology'
    else:
        mapping = {group: group+'-pathology' for group in test_groups}
        adata.obs['disease_group'] = adata.obs['pathology.group'].map(mapping)
        adata.obs.loc[adata.obs['pathology.group'] == reference_group, 'disease_group'] = reference_group+'-pathology'

    # Return the annotated AnnData object
    return adata



def wilcoxon_de(adata: AnnData, 
                reference_group: str = 'no',
                norm_method: str = 'actionet',
                filter: bool = True,
                filter_by: str = 'prop',
                filter_thresh: float = 0.05,
                test_layer: str = 'counts',
                grouped_test: bool = True):
    """
    Perform Wilcoxon rank-sum test for differential expression analysis on cell-type specific subsets of data in AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    norm_method : str, optional
        Normalization method to use. Options are 'actionet', 'scanpy', and 'custom'. Default is 'actionet'.
    filter : bool, optional
        Whether to filter out genes not expressed in a minimum number of cells. Default is True.
    filter_by : str, optional
        Metric to use for filtering. Options are 'ncells' or 'prop'. Default is 'prop'.
    filter_thresh : float, optional
        Threshold for filtering out genes. Either minimum number of cells or proportion of cells expressing gene. Default is 0.05.
    test_layer : str, optional
        Layer to use for differential expression testing. Default is 'counts'.
    reference_group: str, optional.
        Name of the reference group for comparison. Default is 'no'.
    grouped_test: bool, optional
        Whether to group tested levels into one group.

        
    Returns
    -------
    adata_sub : dict
        Dictionary of AnnData objects containing normalized and filtered data for each cell type.

    """

    # Initialize log1p key in adata.uns
    adata.uns['log1p'] = {'base': None}

    # Initialize empty dictionary to store cell-type specific AnnData objects
    adata_sub = {}

    # Iterate over unique cell types in adata and perform differential expression analysis
    for cell_type in adata.obs.cell_type.unique():
        print()

        # Subset data for specific cell type
        adata_sub[cell_type] = adata[adata.obs.cell_type == cell_type].copy()

        # Use counts layer for differential expression testing
        adata_sub[cell_type].X = adata_sub[cell_type].layers['counts'] 

        # Filter out genes not expressed in minimum number of cells or below proportion threshold
        if filter:
            if filter_by=='ncells':
                sc.pp.filter_genes(adata_sub[cell_type], min_cells=10)
            elif filter_by=='prop':
                n_cells = filter_thresh*adata_sub[cell_type].n_obs
                sc.pp.filter_genes(adata_sub[cell_type], min_cells=n_cells)

        # Normalize gene expression data
        if norm_method=='actionet':
            # Normalize using actionet normalization method
            adata_sub[cell_type] = pathway_analyses.normalize_actionet(adata_sub[cell_type], 
                                                                    layer_key = 'counts', layer_key_out = None,
                                                                    top_features_frac = 1.0, scale_factor = "median",
                                                                    transformation = "log", anchor_features = None, copy = True)
        elif norm_method=='scanpy':
            # Normalize using scanpy normalization method
            sc.pp.normalize_total(adata_sub[cell_type])
            sc.pp.log1p(adata_sub[cell_type], layer='counts')

        elif norm_method=='custom':
            # Normalize using custom normalization method
            adata_sub[cell_type] = pathway_analyses.normalize_default(adata_sub[cell_type], log_scale=True)
        
        # Perform Wilcoxon rank-sum test for differential expression analysis        
        print(f'evaluating differential expression in {cell_type}...')
        adata_sub[cell_type].obs['disease_group'] = adata_sub[cell_type].obs['disease_group'].astype('category', copy=False)
        groups = 'all' if grouped_test else [grp for grp in list(adata_sub[cell_type].obs['disease_group'].unique()) if str(grp)!='nan']
        sc.tl.rank_genes_groups(adata_sub[cell_type],
                                method="wilcoxon",
                                use_raw=False,
                                layer=test_layer,
                                groupby = 'disease_group', 
                                reference = reference_group+'-pathology',
                                groups = groups,
                                key_added = 'wilcoxon_test_pathology',
                                corr_method = 'benjamini-hochberg', 
                                verbose = False)
        
    return adata_sub


def get_degs(adata_sub: dict,
             grouped_test: bool = True,
             reference_group: str = 'no',
             test_groups: list = ['late', 'early'],
             save_prefix: str = 'mathys_pfc') -> dict:
    """
    Perform differential gene expression analysis using t-tests and Wilcoxon tests.

    Parameters:
        adata_sub (dict): A dictionary with keys representing cell types and values representing AnnData objects.
        grouped_test (bool): If True, perform Wilcoxon tests between AD and no AD samples for each cell type.
                             If False, perform t-tests between each test group and the reference group for each cell type.
        reference_group (str): The name of the reference group used in t-tests.
        test_groups (list): A list of test group names used in t-tests.
        save_prefix (str): A prefix added to the output file names.

    Returns:
        A dictionary with keys representing test names and cell types, and values representing
        Pandas DataFrames containing the results of differential expression analysis.

    """
    
    degs_t_test = {}
    keys = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']

    if grouped_test:
        degs_t_test['ad_vs_no'] = {}
        with pd.ExcelWriter(f"../results/ad_vs_no/{save_prefix}_t_test_degs.xlsx") as writer:
            for cell_type in adata_sub.keys():
                value_dict = {k:adata_sub[cell_type].uns['wilcoxon_test_pathology'][k].tolist() for k in keys}
                degs_t_test[cell_type] = pd.DataFrame(value_dict)

                degs_t_test[cell_type][['names']] = degs_t_test[cell_type][['names']].applymap(lambda x: x[0])
                degs_t_test[cell_type][['scores', 'pvals', 'pvals_adj', 'logfoldchanges']] = \
                    degs_t_test[cell_type][['scores', 'pvals', 'pvals_adj', 'logfoldchanges']].applymap(lambda x: float(x[0]))

                degs_t_test[cell_type]['abs_logfoldchanges'] = abs(degs_t_test[cell_type]['logfoldchanges'])
                degs_t_test[cell_type].sort_values(by='pvals_adj', inplace=True)

                degs_t_test[cell_type]['Direction_t_test'] = degs_t_test[cell_type]['logfoldchanges'].apply(lambda x:\
                                                    "up" if x>0 else "down")

                degs_t_test[cell_type].reset_index(inplace=True)
                degs_t_test[cell_type].drop('index', axis=1, inplace=True)

                degs_t_test[cell_type].to_excel(writer, sheet_name=cell_type, na_rep='NA')
            # writer.close()
    else:
        for test_name in test_groups:
            temp_test_name = test_name
            test_name = test_name+'_vs_'+reference_group
            degs_t_test[test_name] = {}
            with pd.ExcelWriter(f"../results/{test_name}/{save_prefix}_t_test_degs.xlsx") as writer:
                for cell_type in adata_sub.keys():
                    value_dict = {k:adata_sub[cell_type].uns['wilcoxon_test_pathology'][k][temp_test_name+'-pathology'].tolist() for k in keys}
                    degs_t_test[test_name][cell_type] = pd.DataFrame(value_dict)
                    # degs_t_test[test_name][cell_type][['names']] = degs_t_test[test_name][cell_type][['names']].applymap(lambda x: x[0])
        
                    # degs_t_test[test_name][cell_type][['scores', 'pvals', 'pvals_adj', 'logfoldchanges']] = \
                    #     degs_t_test[test_name][cell_type][['scores', 'pvals', 'pvals_adj', 'logfoldchanges']].applymap(lambda x: float(x[0]))

                    degs_t_test[test_name][cell_type]['abs_logfoldchanges'] = abs(degs_t_test[test_name][cell_type]['logfoldchanges'])
                    degs_t_test[test_name][cell_type].sort_values(by='pvals_adj', inplace=True)

                    degs_t_test[test_name][cell_type]['Direction_t_test'] = degs_t_test[test_name][cell_type]['logfoldchanges'].apply(lambda x:\
                                                        "up" if x>0 else "down")

                    degs_t_test[test_name][cell_type].reset_index(inplace=True)
                    degs_t_test[test_name][cell_type].drop('index', axis=1, inplace=True)

                    degs_t_test[test_name][cell_type].to_excel(writer, sheet_name=cell_type, na_rep='NA')
                # writer.close()

    return degs_t_test         
        