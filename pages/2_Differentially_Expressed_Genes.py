#!/usr/bin/env python
import pandas as pd
import streamlit as st

from constants import DATA
from utils import fetch_resource

st.set_page_config(layout="wide")

st.markdown("""
    # Differentially expressed genes
    
    Genes below have been determined using `sc.tl.rank_genes_groups` where `t-test_overestim_var`
    is the default method.
    """
)

filter_flag = []
ds = st.sidebar.selectbox(
    "**Select models**",
    DATA.keys(),
    index=None,
    placeholder="Select species",
)

if ds:
    filter_by = st.sidebar.selectbox(
        "**Select by**",
        DATA[ds]["DEGS"].keys(),
        index=None,
        placeholder="Select by",
    )

if ds and filter_by:
    markers = pd.read_feather(fetch_resource(DATA[ds]['DEGS'][filter_by]))

    group = st.sidebar.multiselect(
        "**Cell type**", markers.group.unique(), placeholder="Select group ..."
    )

    genes = st.sidebar.multiselect(
        "**Gene**", markers.gene_symbol.unique(), placeholder="Select genes ..."
    )

    foldchange = st.sidebar.number_input(
        "**Log2 fold-change**",
        value=1,
    )

    pval = st.sidebar.number_input(
        "**Adjusted p-value**",
        value=0.05,
    )

    if group:
        filter_flag.append("group == @group")

    if genes:
        filter_flag.append("@genes in gene_symbol")

    if foldchange:
        filter_flag.append(
            "logfoldchanges > @foldchange"
            if foldchange > 0
            else "logfoldchanges < @foldchange"
        )

    if pval:
        filter_flag.append("pvals_adj < @pval")

    subset = markers.query(" & ".join(filter_flag)) if filter_flag else markers
    st.dataframe(subset, use_container_width=True, height=650)
