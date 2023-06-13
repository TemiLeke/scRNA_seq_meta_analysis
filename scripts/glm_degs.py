import warnings
from collections import Counter

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import bs, cr, dmatrix
from scipy import stats
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
from tqdm import tqdm




def glm_degs(
    adata,
    X_data=None,
    genes=None,
    layer=None,
    fullModelFormulaStr="~disease_group + n_genes + Subjects",
    reducedModelFormulaStr="~1",
    family="NegativeBinomial",
    ):

    """Differential genes expression tests using generalized linear regressions.

    Tests each gene for differential expression as a function of specified covariates (e.g the disease group/condition in `adata.ob`)
    using generalized additive models. This function can
    also use other covariates as specified in the full (i.e `~clusters`) and reduced model formula to identify differentially
    expression genes across different categories, group, etc.

    glm_degs relies on statsmodels package and is adapted from the `differentialGeneTest` function in Monocle. Note that
    glm_degs supports performing deg analysis for any layer or normalized data in your adata object.


    Parameters
    ----------
        adata: :class:`~anndata.AnnData`
            an Annodata object
        X_data: `np.ndarray` (default: `None`)
            The user supplied data that will be used for differential expression analysis directly.
        genes: `list` or None (default: `None`)
            The list of genes that will be used to subset the data for differential expression analysis. If `None`, all
            genes will be used.
        layer: `str` or None (default: `None`)
            The layer that will be used to retrieve data for dimension reduction and clustering. If `None`, .X is used.
        fullModelFormulaStr: `str` (default: `~disease_group + n_genes + Subjects`)
            A formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for
            each gene/feature.
        reducedModelFormulaStr: `str` (default: `~1`)
            A formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests)
            for each gene/feature.
        family: `str` (default: `NB2`)
            The distribution family used for the expression responses in statsmodels. Currently always uses `NB2` and this
            is ignored. NB model requires us to define a parameter $\alpha$ which it uses to express the variance in terms
            of the mean as follows: variance = mean + $\alpha$ mean^p. When $p=2$, it corresponds to the NB2 model. In order
            to obtain the correct parameter $\alpha$ (sm.genmod.families.family.NegativeBinomial(link=None, alpha=1.0), by
            default it is 1), we use the auxiliary OLS regression without a constant from Messrs Cameron and Trivedi. More
            details can be found here: https://towardsdatascience.com/negative-binomial-regression-f99031bb25b4.

    Returns
    -------
        Returns an updated `~anndata.AnnData` with a new key `glm_degs` in the .uns attribute, storing the differential
        expression test results after the GLM test.
    """

    if X_data is None:
        if layer is None:
            X_data = adata.X
        else:
            X_data = adata.layers[layer]

    if genes is None:
        genes = list(adata.var_names)

    assert len(genes) == X_data.shape[1], f"When providing X_data, a \
                                            list of genes name that corresponds \
                                            to the columns of X_data must be provided"

                
    if layer is None:
        if issparse(X_data):
            X_data.data = 2**X_data.data - 1
        else:
            X_data = 2**X_data - 1

    factors = get_all_variables(fullModelFormulaStr)
    factors = [i for i in factors]

    if len(set(factors).difference(adata.obs.columns)) == 0:
        df_factors = adata.obs[factors].copy()
    else:
        raise Exception(
            f"adata object doesn't include the factors from the model formula " f"{fullModelFormulaStr} you provided."
        )

    sparse = issparse(X_data)
    deg_df = pd.DataFrame(index=genes, columns=["status", "family", "pval"])
    
    df_factors['expression'] = None

    for i, gene in tqdm(enumerate(genes), "Detecting disease category dependent genes via Generalized Additive Models (GAMs)",):
        expression = X_data[:, i].A.copy() if sparse else X_data[:, i].copy()
        df_factors["expression"] = expression
        deg_df.iloc[i, :] = diff_test_helper(df_factors, fullModelFormulaStr, reducedModelFormulaStr, family)
        df_factors['expression'] = None

    deg_df["qval"] = multipletests(deg_df["pval"], method="fdr_bh")[1]

    adata.uns["glm_degs"] = deg_df



def diff_test_helper(
    data,
    fullModelFormulaStr="~disease_group + n_genes",
    reducedModelFormulaStr="~1",
    family='NegativeBinomial',
):
    # Dividing data into train and validation datasets
    transformed_x = dmatrix(fullModelFormulaStr, data, return_type="dataframe")
    transformed_x_null = dmatrix(reducedModelFormulaStr, data, return_type="dataframe")

    expression = data["expression"]

    try:
        # poisson_training_results = sm.GLM(expression, transformed_x, family=sm.families.Poisson()).fit()
        # poisson_df = pd.DataFrame({"mu": poisson_training_results.mu, "expression": expression})
        # poisson_df["AUX_OLS_DEP"] = poisson_df.apply(
        #     lambda x: ((x["expression"] - x["mu"]) ** 2 - x["expression"]) / x["mu"],
        #     axis=1,
        # )
        # ols_expr = """AUX_OLS_DEP ~ mu - 1"""
        # aux_olsr_results = smf.ols(ols_expr, poisson_df).fit()

        nb2_family = eval(f'sm.families.{family}()')  # (alpha=aux_olsr_results.params[0])

        nb2_full = sm.GLM(expression, transformed_x, family=nb2_family).fit()
        nb2_null = sm.GLM(expression, transformed_x_null, family=nb2_family).fit()
    except:
        return ("fail", "NB2", 1)

    pval = lrt(nb2_full, nb2_null)
    
    return ("ok", "NB2", pval)


def get_all_variables(formula):
    md = patsy.ModelDesc.from_formula(formula)
    termlist = md.rhs_termlist + md.lhs_termlist

    factors = []
    for term in termlist:
        for factor in term.factors:
            factors.append(factor.name())

    return factors


def lrt(full, restr):
    llf_full = full.llf
    llf_restr = restr.llf
    df_full = full.df_resid
    df_restr = restr.df_resid
    lrdf = df_restr - df_full
    lrstat = -2 * (llf_restr - llf_full)
    lr_pvalue = stats.chi2.sf(lrstat, df=lrdf)

    return lr_pvalue





    # adata_exc.obs.amyloid_x = adata_exc.obs.amyloid_x.astype(float)
# adata_exc.obs.plaq_n_x = adata_exc.obs.plaq_n_x.astype(float)
# adata_exc.obs.age_death.cat.rename_categories({"90+": '90'}, inplace=True)
# adata_exc.obs.age_death = adata_exc.obs.age_death.astype(float)
# adata_exc.obs.nft = adata_exc.obs.nft.astype(float)
# adata_exc.obs.tangles = adata_exc.obs.tangles.astype(float)

# adata_exc.obs['ADcond'] = adata_exc.obs['pathologic diagnosis of AD'].copy()


%%R -i adata_exc -o exc_degs

##############################################################################################################
library(nebula)


adata_exc$amyloid <- adata_exc$amyloid_x
adata_exc$plaq <- adata_exc$plaq_n_x
adata_exc$gpath <- adata_exc$gpath_x
adata_exc$sex <- adata_exc$msex_x

##############################################################################################################
Predictors <- as.data.frame(adata_exc@colData[,c("amyloid", "plaq", "sex", "age_death", "ADcond", 'n_genes')])

print('scaling numerical factors...')
Predictors$amyloid <- scale(Predictors$amyloid)
Predictors$age_death <- scale(Predictors$age_death)
Predictors$plaq <- scale(Predictors$plaq)

print('converting categorical factors...')
Predictors$sex <- as.character(Predictors$sex)
Predictors$ADcond <- as.character(Predictors$ADcond)

# build design matrix from a
df = model.matrix(~ADcond + amyloid + n_genes, data=Predictors)

print('grouping cells...')
# group cells of the same Subject/patient id

data_g <- group_cell(count=assays(adata_exc)$counts, id=adata_exc$Subject, pred=df)
sparsematrix <- as(data_g$count, "CsparseMatrix")


print('computing degs...')
re = nebula(sparsematrix, data_g$id, pred=data_g$pred, offset=Matrix::colSums(sparsematrix), method='HL', model='PMM')
##############################################################################################################
exc_degs <- re$summary
rownames(exc_degs) <- exc_degs$gene


exc_degs$score <- -log10(exc_degs$'p_ADcondYES') * ifelse(exc_degs$'p_ADcondYES'>0, 1, -1)
print('saving degs...')

exc_degs = exc_degs[,c('p_ADcondYES', 'logFC_ADcondYES', 'score','gene')]

gc()

colnames(exc_degs) = c('p_val', 'logFC_AD_vs_NO', 'score', 'gene')

# saveRDS(out, file="../data/differentially_expressed_genes/E4_nebula_associations_by_celltype_Oli.rds")
##############################################################################################################
