#!/usr/bin/env python
import anndata
import streamlit as st

from constants import DATA, DEFAULT_DR, DEFAULT_META
from utils import fetch_resource, plot_feature, plot_sc_embedding

st.set_page_config(layout="wide")
st.markdown("""
    # Gene expression
    
    Levels of gene activity along differentiation.
""")

ds = st.sidebar.selectbox(
    "**Load dataset**",
    DATA.keys(),
    index=None,
    placeholder="Select contact method...",
)

if ds is not None:

    adata = anndata.read_h5ad(fetch_resource(DATA[ds]['DATASET']))

    sl_dr = st.sidebar.selectbox(
        "**Dimension reduction**",
        adata.obsm_keys(),
        index=adata.obsm_keys().index(DEFAULT_DR),
        placeholder="Select method ...",
    )

    sl_metadata = st.sidebar.selectbox(
        "**Metadata**",
        adata.obs.columns,
        index=adata.obs.columns.get_loc(DEFAULT_META),
        placeholder="Select column ...",
    )

    sl_feature = st.sidebar.selectbox(
        "**Gene**",
        adata.raw.var_names,
        index=0,
        placeholder="Select gene ...",
    )

    is_imputed = sl_feature in adata.var_names
    sl_denoised = st.sidebar.checkbox(
        "Use denoised expression?",
        help="Denoised expression is sampled from the decoder.",
        disabled=(not is_imputed)
    )

    col1, col2 = st.columns(2)
    plot_sc_embedding(adata, group_by=sl_metadata, reduction_key=sl_dr, ax=col1)
    plot_sc_embedding(
        adata, feature=sl_feature, reduction_key=sl_dr, layer=sl_denoised, ax=col2
    )

    st.markdown("## Raw gene expression")
    plot_feature(adata, feature=sl_feature, group_by=sl_metadata, kind="box")
