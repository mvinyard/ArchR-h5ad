
__module_name__ = "_cleanup_anndata.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import scipy.sparse


format_dict = {
    "csc": scipy.sparse.csc_matrix,
    "csr": scipy.sparse.csr_matrix,
}

def _to_sparse_format(X, to="csr"):
    return format_dict[to](X)

def _format_adata_indices(adata):
    
    adata.obs.index = adata.obs.index.astype(str)
    adata.var.index = adata.var.index.astype(str)
    
    return adata

def _cleanup_anndata(adata):
    
    adata.X = _to_sparse_format(adata.X, to="csr")
    return _format_adata_indices(adata)