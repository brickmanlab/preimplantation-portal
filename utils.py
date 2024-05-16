import anndata
import pandas as pd
import plotly.express as px
import streamlit as st

VERSION = 1.1
ZENODO_URL = "https://zenodo.org/records/11204495/files"
DATA = {
    "human": {
        "ds": f"{ZENODO_URL}/portal_human_v1.1.h5ad",
        "degs_ct": f"{ZENODO_URL}/human_degs_ct.feather",
        "degs_stage": f"{ZENODO_URL}/human_degs_stage.feather",
        "shap": f"{ZENODO_URL}/human_SHAP.feather",
    },
    "mouse": {
        "ds": f"{ZENODO_URL}/portal_mouse_v1.1.h5ad",
        "degs_ct": f"{ZENODO_URL}/mouse_degs_ct.feather",
        "degs_stage": f"{ZENODO_URL}/mouse_degs_stage.feather",
        "shap": f"{ZENODO_URL}/mouse_SHAP.feather",
    },
}


def get_embedding(adata: anndata.AnnData, key: str) -> pd.DataFrame:
    """
    Helper function which retrieves embedding coordinates for each cell.

    Parameters
    ----------
    adata : anndata.AnnData
        scrna-seq dataset
    key : str
        Dimension reduction key, usually starts with X_

    Returns
    -------
    pd.DataFrame
        Embedding coordinates

    Raises
    ------
    ValueError
        Fail if reduction key doesn't exist
    """
    if key not in adata.obsm.keys():
        raise ValueError(f"Reduction key: {key} not available")

    dimension_names = f"{key[2:].upper()}_1", f"{key[2:].upper()}_2"
    return pd.DataFrame(adata.obsm[key][:, :2], columns=dimension_names)


def plot_sc_embedding(
    adata: anndata.AnnData,
    reduction_key: str,
    group_by: str = None,
    feature: str = None,
    layer: str = None,
    ax=None,
):
    """
    Plot single-cell dataset

    Parameters
    ----------
    adata : anndata.AnnData
        scrna-seq dataset
    reduction_key : str
        Reduced space key
    group_by : str
        Key used to color cells
    features: str
        Gene
    ax : _type_
        Axes
    """
    embeddings = get_embedding(adata, reduction_key)

    if group_by:
        embeddings[group_by] = adata.obs[group_by].values
        embeddings = embeddings.sort_values(by=group_by)

        # color_uns_key = f"{group_by}_colors"

        kwargs = {"color": embeddings[group_by].values.tolist()}
        if adata.obs[group_by].dtype == "category":
            ...
        else:
            kwargs["color_continuous_scale"] = px.colors.sequential.Viridis

    if feature:
        X = (
            adata[:, feature].layers["scVI_normalized"].A
            if layer
            else adata.raw[:, feature].X.A
        )
        embeddings[feature] = X.ravel()
        kwargs = {
            "color": embeddings[feature].values.tolist(),
            # "title": feature,
            "color_continuous_scale": px.colors.sequential.Viridis,
        }

    ax_ = ax if ax else st
    ax_.plotly_chart(
        px.scatter(
            data_frame=embeddings,
            x=embeddings.columns[0],
            y=embeddings.columns[1],
            **kwargs,
        ),
        use_container_width=True,
        # .update_xaxes(showgrid=False)
        # .update_yaxes(showgrid=False, zeroline=False)
    )
