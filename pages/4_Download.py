#!/usr/bin/env python
import streamlit as st

st.set_page_config(layout="wide")

st.markdown(
    """
    # Download

    ## Pipelines

    - https://github.com/nf-core/fetchngs/ (revision 1.10.0)
    - https://github.com/brickmanlab/scrnaseq (revision: feature/smartseq)
    

    ## Data analysis
    
    - https://github.com/brickmanlab/preimplantation-models

    Raw pre-processed files

    - Zenodo Link #1
    - Zenodo Link #2

    ## AI models

    Trained models with parameters were uploaded to Hugging Face
    
    - https://huggingface.co/brickmanlab/preimplantation-models

    ### Models

    - scANVI
    - scANVI [ns=15]
    - XGBoost
        - scVI
        - scANVI
        - scGEN

    """
)
