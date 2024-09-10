#!/usr/bin/env python
import streamlit as st

st.set_page_config(layout="wide")

st.markdown(
    """
    # Deep Learning Based Models for Preimplantation Mouse and Human Embryos Based on Single Cell RNA Sequencing

    _[Martin Proks](https://orcid.org/0000-0002-8178-3128)\*, 
    [Nazmus Salehin](https://orcid.org/0000-0002-8155-4296)\*, 
    [Joshua M. Brickman](https://orcid.org/0000-0003-1580-7491)**_
    
    _\* There authors contributed equally to the work_

    _** Corresponding author [joshua.brickman@sund.ku.dk](mailto:joshua.brickman@sund.ku.dk)_

    The rapid growth of single-cell transcriptomic technology has produced an increasing number of
    datasets for both embryonic development and _in vitro_ pluripotent stem cell derived models.
    This avalanche of data surrounding pluripotency and the process of lineage specification has
    meant it has become increasingly difficult to define specific cell types or states and compare
    these to _in vitro_ differentiation. Here we utilize a set of deep learning (DL) tools to
    integrate and classify multiple datasets. This allows for the definition of both mouse and
    human embryo cell types, lineages and states, thereby maximising the information one can garner
    from these precious experimental resources. Our approaches are built on recent initiatives for
    large scale human organ atlases, but here we focus on the difficult to obtain and process
    material that spans early mouse and human development. We deploy similar approaches as the
    initiatives building large reference organ atlases, however with a focus on early mammalian
    development. Using publicly available data for these stages, we test different deep learning
    approaches and develop a model to classify cell types in an unbiased fashion at the same time as
    defining the set of genes used by the model to identify lineages, cell types and states. We have
    used our models trained on _in vivo_ development to classify pluripotent stem cell models for
    both mouse and human development, showcasing the importance of this resource as a dynamic
    reference for early embryogenesis.
    """
)

st.image(
    "static/Fig-1.v4.3.png",
    caption="""
        Summary of datasets used to build reference models. a) Schematic overview of mouse
        and human preimplantation development. b) Quantification of cells per publication which
        were collected for building the mouse (grey) and human (black) reference. c) Computational
        schematic of tools used to build and interpret the reference models. d) Gene expression
        of canonical markers for each developmental stage in mouse (top) and human (bottom)
        preimplantation. e) Reduced dimensional representation of preimplantation mouse (left)
        and human (right) datasets. dpf: days post fertilization, E: embryonic day.
         """,
)
