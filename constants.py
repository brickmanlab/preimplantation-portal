VERSION = 1.5

DEFAULT_DR = "X_draw_graph_fa"
DEFAULT_META = "stage"
ZENODO_URL = "https://zenodo.org/records/13749348/files"

DATA = {
    "HUMAN": {
        "RAW_DATASET": f"{ZENODO_URL}/01_mouse_reprocessed.h5ad",
        "DATASET": f"{ZENODO_URL}/portal_human_v{VERSION}.h5ad",
        "DEGS": {
            "CT": f"{ZENODO_URL}/human_degs_ct_v{VERSION}.feather",
            "STAGE": f"{ZENODO_URL}/human_degs_stage_v{VERSION}.feather"
        },
        "SHAP": f"{ZENODO_URL}/human_SHAP_v{VERSION}.feather",
    },
    "MOUSE": {
        "RAW_DATASET": f"{ZENODO_URL}/32_human_adata.h5ad",
        "DATASET": f"{ZENODO_URL}/portal_mouse_v{VERSION}.h5ad",
        "DEGS": {
            "CT": f"{ZENODO_URL}/mouse_degs_ct_v{VERSION}.feather",
            "STAGE": f"{ZENODO_URL}/mouse_degs_stage_v{VERSION}.feather",
        },
        "SHAP": f"{ZENODO_URL}/mouse_SHAP_v{VERSION}.feather",
    },
}
