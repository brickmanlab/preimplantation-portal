#!/usr/bin/env python
import streamlit as st

st.set_page_config(layout="wide")

st.markdown(
    """
    # Deep Learning Based Models for Preimplantation Mouse and Human Development

    _[Martin Proks](https://orcid.org/0000-0002-8178-3128)\*, 
    [Nazmus Salehin](https://orcid.org/0000-0002-8155-4296)\*, 
    [Joshua M. Brickman](https://orcid.org/0000-0003-1580-7491)**_
    
    _\* There authors contributed equally to the work_
    
    _** Corresponding author [joshua.brickman@sund.ku.dk](mailto:joshua.brickman@sund.ku.dk)_

    Abstract goes here!
    """
)

st.image(
    "static/Fig-1.v4.3.png",
    caption="""
         Summary of datasets used to build reference models. A) Schematic overview of mouse and human 
         preimplantation development. B) Quantification of cells per publication which were collected 
         for building the mouse (grey) and human (black) reference. C) Computational schematic of 
         tools used to build and interpret the reference models. D) Reduced dimensional representation 
         of preimplantation mouse (left) and human (right) datasets. E) Gene expression of canonical 
         markers for each developmental stage in mouse (top) and human (bottom) preimplantation.
         """,
)
