#!/usr/bin/env python
import urllib.request
from pathlib import Path

import anndata
import pandas as pd
import streamlit as st

from utils import DATA, plot_sc_embedding

st.set_page_config(layout="wide")
DEFAULT_DR = "X_draw_graph_fa"
DEFAULT_META = "stage"


@st.cache_resource
def load_dataset(ds: str) -> anndata.AnnData:
    default_path: str = f"/tmp/{ds}.h5ad"
    if not Path(default_path).exists():
        with st.spinner("Please wait we are downloading the SDH Model."):
            urllib.request.urlretrieve(DATA[ds]["ds"], default_path)

    return anndata.read_h5ad(default_path)


def plot_embedding(
    adata: anndata.AnnData, dr_key: str, groupby: str, dims: int = 2, ax=None
):
    if dr_key not in adata.obsm_keys():
        st.error(f"{dr_key} not present in dataset", icon="🚨")

    df = pd.DataFrame(
        adata.obsm[dr_key][:, :dims],
        columns=[f"{dr_key}_{i}" for i in range(1, dims + 1)],
    )
    df[groupby] = adata.obs[groupby].values

    color_scheme = (
        adata.uns[f"{groupby}_colors"]
        if adata.obs[groupby].dtype.name == "category"
        else "viridis"
    )

    print(color_scheme)

    # selection = alt.selection_point(fields=[groupby], bind="legend")
    # chart = (
    #     (
    #         alt.Chart(df)
    #         .mark_circle()
    #         .encode(
    #             x=df.columns[0],
    #             y=df.columns[1],
    #             color=groupby,
    #             tooltip=[groupby],
    #             opacity=alt.condition(selection, alt.value(1), alt.value(0.05)),
    #         )
    #         .interactive()
    #     )
    #     # .configure_range(category=color_scheme)
    #     .configure_range(category={"scheme": "Dark2"})
    #     .configure_axis(grid=False, labels=False)
    #     .add_params(selection)
    # )

    # ax_ = ax if ax is not None else st
    # ax_.altair_chart(chart, theme=None, use_container_width=True)


st.markdown("# Gene expression")

ds = st.sidebar.selectbox(
    "**Load dataset**",
    DATA.keys(),
    index=None,
    placeholder="Select contact method...",
)

if ds is not None:
    adata = load_dataset(ds)

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
        "Use denoised expression?", disabled=(not is_imputed)
    )

    col1, col2 = st.columns(2)
    plot_sc_embedding(adata, group_by=sl_metadata, reduction_key=sl_dr, ax=col1)
    plot_sc_embedding(
        adata, feature=sl_feature, reduction_key=sl_dr, layer=sl_denoised, ax=col2
    )

    # st.markdown("## Heatmap")

    # sl_features = st.selectbox(
    #     "**Select multiple genes**", adata.var_names, placeholder="Select genes ..."
    # )

    # if msl_metadata is not None and sl_dr is not None:
    #     tabs = st.columns(len(msl_metadata))
    #     for idx, m in enumerate(msl_metadata):
    #         plot_embedding(adata, dr_key=sl_dr, groupby=m, ax=tabs[idx])
