#!/usr/bin/env python
import streamlit as st

st.set_page_config(layout="wide")

st.markdown(
    """
    # Download

    ## 1. Pipelines

    - Downloading datasets: [nf-core/fetchngs (revision 1.10.0)](https://github.com/nf-core/fetchngs)
    - Aligning datasets: [brickmanlab/scrnaseq (revision: feature/smartseq)](https://github.com/brickmanlab/scrnaseq)
    - **Ensembl Genomes**
        - Mouse: GRCm38 v102
        - Human: GRCh38 v110

    ## 2. Codebase
    
    - GitHub: [brickmanlab/preimplantation-models](https://github.com/brickmanlab/preimplantation-models)
    - Portal codebase: [brickmanlab/preimplantation-portal](https://github.com/brickmanlab/preimplantation-portal)
    
    ## 3. Raw data
    
    - [Mouse](https://zenodo.org/records/11204495/files/01_mouse_reprocessed.h5ad?download=1)
    - [Human](https://zenodo.org/records/11204495/files/32_human_adata.h5ad?download=1)


    ## 4. AI models

    Trained models with parameters were uploaded to [Hugging Face](https://huggingface.co/brickmanlab/preimplantation-models).

    ### 4.1 Models

    - [scANVI mouse](https://huggingface.co/brickmanlab/mouse-scanvi)
    - [scANVI human](https://huggingface.co/brickmanlab/human-scanvi)

    """
)
