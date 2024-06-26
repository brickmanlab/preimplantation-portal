#!/usr/bin/env python
import pandas as pd
import streamlit as st

from utils import DATA

st.set_page_config(layout="wide")

st.markdown("# SHAP features")

ds = st.sidebar.selectbox(
    "**Load dataset**",
    DATA.keys(),
    index=None,
    placeholder="Select dataset ...",
)

if ds:
    data = pd.read_feather(DATA[ds]["shap"])
    subset = data.copy()

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
        subset = subset.query(" & ".join(filter_condition))

    st.dataframe(subset, use_container_width=True, height=650)
