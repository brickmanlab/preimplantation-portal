#!/usr/bin/env python
import streamlit as st

from constants import DATA

st.set_page_config(layout="wide")

st.markdown(
    f"""
    # Download

    ## 1. Pipelines

    - Downloading datasets: [nf-core/fetchngs (revision 1.10.0)](https://github.com/nf-core/fetchngs)
    - Aligning datasets: [brickmanlab/scrnaseq (revision: feature/smartseq)](https://github.com/brickmanlab/scrnaseq)
    - **Ensembl Genomes**
        - Mouse: GRCm38 v102
        - Human: GRCh38 v110

    ## 2. Codebase
    
    - Data analysis: [brickmanlab/proks-salehin-et-al](https://github.com/brickmanlab/proks-salehin-et-al)
    - Web portal: [brickmanlab/preimplantation-portal](https://github.com/brickmanlab/preimplantation-portal)
    
    ## 3. Raw data
    
    - [Mouse]({DATA['MOUSE']['RAW_DATASET']})
    - [Human]({DATA['HUMAN']['RAW_DATASET']})

    ## 4. AI models

    Trained models with parameters were uploaded to [Hugging Face](https://huggingface.co/brickmanlab/preimplantation-models).

    ### 4.1 Models

    - [scANVI mouse](https://huggingface.co/brickmanlab/mouse-scanvi)
    - [scANVI human](https://huggingface.co/brickmanlab/human-scanvi)

    """
)
