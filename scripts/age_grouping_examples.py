# scripts/age_grouping_examples.py
# Usage: edit DATA_FILE to point to your CSV; then run with python scripts/age_grouping_examples.py
# Expects CSV with rows = samples and columns: proteins named starting with 'p', and columns: age, sex, batch, CRP (optional), SampleID (optional)

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.api import OLS
from statsmodels.tools.tools import add_constant
from statsmodels.stats.multitest import multipletests

SEED = 12345
np.random.seed(SEED)

# === User params ===
DATA_FILE = 'data/expr_sample_table.csv'  # change to your SomaScan/Olink CSV path
OUT_DIR = 'outputs'
PROT_PREFIX = 'p'   # set to '' to auto-detect protein columns (recommended for real SomaScan/Olink)
PROTEIN_LIST_FILE = 'data/top50_proteins.txt'  # optional, one protein per line; if present, used directly
TOP_N = 50  # number of top proteins to select for heatmap
REF_LEVEL = 'Q1'  # reference for quartiles
# ====================

os.makedirs(OUT_DIR, exist_ok=True)

print('Loading data:', DATA_FILE)
df = pd.read_csv(DATA_FILE)

# quartiles and custom bins
q_breaks = df['age'].quantile([0,0.25,0.5,0.75,1.0]).values
# avoid duplicate breaks (constant age) by jittering tiny amount if needed
if len(np.unique(q_breaks)) < len(q_breaks):
    q_breaks = np.linspace(df['age'].min()-1e-6, df['age'].max()+1e-6, 5)

df['age_q'] = pd.qcut(df['age'], q=4, labels=['Q1','Q2','Q3','Q4'])
# custom example
bins = [-np.inf, 50, 65, 80, np.inf]
labels = ['<50','50-64','65-79','80+']
df['age_grp'] = pd.cut(df['age'], bins=bins, labels=labels, include_lowest=True)
# ordinal coding
if pd.api.types.is_categorical_dtype(df['age_q']):
    df['age_ord'] = df['age_q'].cat.codes + 1
else:
    df['age_ord'] = pd.Categorical(df['age_q']).codes + 1

# protein columns
if PROT_PREFIX == '':
    # auto-detect protein columns: numeric columns excluding known metadata
    metadata_cols = set(['age','sex','batch','CRP','SampleID','age_q','age_grp','age_ord'])
    prot_cols = [c for c in df.columns if c not in metadata_cols and pd.api.types.is_numeric_dtype(df[c])]
else:
    prot_cols = [c for c in df.columns if c.startswith(PROT_PREFIX)]
if len(prot_cols) == 0:
    raise SystemExit('No protein columns found; set PROT_PREFIX or provide a protein list file')

# covariates
covar_cols = []
for c in ['sex','batch','CRP']:
    if c in df.columns:
        covar_cols.append(c)

# prepare base design matrix (one-hot encode categorical covariates drop_first)
X_base = pd.DataFrame(index=df.index)
for c in covar_cols:
    # treat non-numeric as categorical
    if not pd.api.types.is_numeric_dtype(df[c]):
        dummies = pd.get_dummies(df[c].astype(str), prefix=c, drop_first=True)
        X_base = pd.concat([X_base, dummies], axis=1)
    else:
        # coerce numeric, allow NaN
        X_base[c] = pd.to_numeric(df[c], errors='coerce')
# add intercept
X_base = add_constant(X_base, has_constant='add')

# age dummies for quartiles (drop_first -> Q1 reference)
age_dummies = pd.get_dummies(df['age_q'].astype(str), prefix='age_q', drop_first=True)
if age_dummies.shape[1] == 0:
    raise SystemExit('age_q dummy not created; check age data')

X_full = pd.concat([X_base, age_dummies], axis=1)

results = []
import statsmodels.api as sm
for prot in prot_cols:
    y = df[prot].astype(float)
    mask = (~y.isna()) & X_full.notnull().all(axis=1)
    if mask.sum() < 3:
        results.append((prot, np.nan, np.nan, np.nan, np.nan, mask.sum()))
        continue
    X = X_full.loc[mask, :]
    ysub = y.loc[mask]
    try:
        model = sm.OLS(ysub.values, X.values).fit()
        # coefficient for highest quartile vs Q1: find age_q_Q4 column
        age_col = [c for c in X.columns if 'age_q_Q4' in c]
        if len(age_col) == 0:
            # fallback to first age_q column
            age_col = [c for c in X.columns if c.startswith('age_q_')]
        if len(age_col) == 0:
            coef = se = t = pval = np.nan
        else:
            idx = list(X.columns).index(age_col[0])
            coef = model.params[idx]
            se = model.bse[idx]
            t = model.tvalues[idx]
            pval = model.pvalues[idx]
    except Exception as e:
        coef = se = t = pval = np.nan
    results.append((prot, coef, se, t, pval, mask.sum()))

res_df = pd.DataFrame(results, columns=['protein','beta','se','t','pvalue','n_nonmiss'])
res_df['adj_p'] = multipletests(res_df['pvalue'].fillna(1), method='fdr_bh')[1]
res_df = res_df.sort_values('adj_p')
res_df.to_csv(os.path.join(OUT_DIR, 'perprotein_ols_results.csv'), index=False)

# select top proteins
# if a protein list file exists, use it (one protein per line)
if os.path.exists(PROTEIN_LIST_FILE):
    with open(PROTEIN_LIST_FILE, 'r', encoding='utf-8') as fh:
        listed = [ln.strip() for ln in fh if ln.strip()]
    sel = [p for p in listed if p in df.columns]
    if len(sel) == 0:
        raise SystemExit('Protein list provided but none found in data columns')
else:
    sel = res_df.dropna(subset=['adj_p']).head(TOP_N)['protein'].values
    sel = [p for p in sel if p in df.columns]

if len(sel) == 0:
    raise SystemExit('No selected proteins found for heatmap')

expr_top = df[sel].copy()
# z-score per protein
expr_top_z = expr_top.apply(lambda col: (col - np.nanmean(col)) / np.nanstd(col, ddof=0), axis=0)

# column colors
lut = {'Q1':'#1f77b4','Q2':'#2ca02c','Q3':'#9467bd','Q4':'#ff7f0e','old':'#ff7f0e','young':'#1f77b4'}
col_colors = df['age_q'].map(lut)

# clustermap (proteins rows)
# use cividis colormap (good for colorblindness, perceptually uniform)
cmap_name = 'cividis'
g = sns.clustermap(expr_top_z.T, cmap=cmap_name, row_cluster=True, col_cluster=True, col_colors=col_colors.values, figsize=(10,8), xticklabels=False)
plt.suptitle('Top {} age-differential proteins (OLS adj covariates)'.format(len(sel)), y=1.02)
plt.savefig(os.path.join(OUT_DIR, 'age_diff_heatmap_ols.png'), dpi=300, bbox_inches='tight')
plt.close()

print('Done. Results and heatmap saved to', OUT_DIR)
