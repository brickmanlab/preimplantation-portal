#!/usr/bin/env python
import pandas as pd
import streamlit as st

from constants import DATA
from utils import fetch_resource

st.set_page_config(layout="wide")

st.markdown("""
    # SHAP features
    
    Predicted features (genes) used by the scANVI classifier to determine a cell type. The features
    have been determined using [SHAP](https://shap.readthedocs.io/en/latest/).
    
    Each metric for a feature is determined from 10 random boostraps with replacement.
            
    - weight_mean: $\mu$ of SHAP value
    - weight_std: $\sigma$ of SHAP value
    - weight_ci_upper: $\mu$ + $\sigma$
    - weight_ci_lower: $\mu$ - $\sigma$
    - logfoldchanges: Log2fold change from differentiation expression analysis
    - pvals_adj: Adjusted p-value from differentiation expression analysis
    - scores: Estimated score from differentiation expression analysis
    """
)

ds = st.sidebar.selectbox(
    "**Load dataset**",
    DATA.keys(),
    index=None,
    placeholder="Select dataset ...",
)

if ds:
    data = pd.read_feather(fetch_resource(DATA[ds]["SHAP"]))

    query = st.sidebar.selectbox(
        "**Subset**",
        data.ct.unique().tolist(),
        index=None,
        placeholder="Select cell type ...",
    )

    features = st.sidebar.multiselect(
        "**Genes**", data.feature.unique(), placeholder="Select genes ..."
    )

    filter_condition = []
    if query:
        filter_condition.append("ct == @query")
    if features:
        filter_condition.append("feature in @features")

    if filter_condition:
        data = data.query(" & ".join(filter_condition))

    st.dataframe(data, use_container_width=True, height=650)
