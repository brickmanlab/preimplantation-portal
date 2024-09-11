VERSION = 1.5

DEFAULT_DR = "X_draw_graph_fa"
DEFAULT_META = "stage"
ZENODO_URL = "https://zenodo.org/records/13749348/files"

DATA = {
    "HUMAN": {
        "DATASET": f"{ZENODO_URL}/portal_human_v1.1.h5ad",
        "DEGS": {
            "CT": f"{ZENODO_URL}/human_degs_ct.feather",
            "STAGE": f"{ZENODO_URL}/human_degs_stage.feather"
        },
        "SHAP": f"{ZENODO_URL}/human_SHAP.feather",
    },
    "MOUSE": {
        "DATASET": f"{ZENODO_URL}/portal_mouse_v1.1.h5ad",
        "DEGS": {
            "CT": f"{ZENODO_URL}/mouse_degs_ct.feather",
            "STAGE": f"{ZENODO_URL}/mouse_degs_stage.feather",
        },
        "SHAP": f"{ZENODO_URL}/mouse_SHAP.feather",
    },
}
