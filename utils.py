import socket
import urllib.request
from pathlib import Path
from typing import Literal

import anndata
import pandas as pd
import plotly.express as px
import streamlit as st


@st.cache_data
def fetch_resource(url: str) -> str:
    """Helper function for downloading datasets

    Parameters
    ----------
    url : str
        Zenodo url link

    Returns
    -------
    str
        Path where the file was downloaded to, default /tmp
    """
    
    filename = f"/tmp/{url.split('/')[-1]}"
    if not Path(filename).exists():
        try:
            urllib.request.urlretrieve(url, filename)
        except (socket.gaierror, urllib.error.URLError) as err:
            raise ConnectionError(f"could not download {url} due to {err}")

    return filename


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
    ax = None,
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
            adata[:, feature].layers["scVI_normalized"].toarray()
            if layer
            else adata.raw[:, feature].X.toarray()
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


def plot_feature(
    adata: anndata.AnnData, 
    feature: str, 
    group_by: str, 
    kind: Literal["box"] = "box", 
    ax = None
):
    """Plot feature expression

    Parameters
    ----------
    adata : anndata.AnnData
        Dataset
    feature : str
        Gene name
    group_by : str
        Metadata column
    kind : str
        Type of plot
    ax : _type_, optional
        Axis, by default None
    """

    df = pd.DataFrame(adata.raw[:, feature].X.toarray(), columns=[feature])
    df[group_by] = adata.obs[group_by].values
    df = df.sort_values(by=group_by)

    g = None
    match(kind):
        case "box":
            g = px.box(df, x=group_by, y=feature, color=group_by)
        case _:
            raise ValueError(f"Provided kind: {kind} not supported")

    ax_ = ax if ax else st
    ax_.plotly_chart(g, use_container_width=True)
