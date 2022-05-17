
__module_name__ = "_compose_anndata.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])




# import local dependencies #
# ------------------------- #
from ._dict_to_adata import _dict_to_adata
from ._add_obs_var import _add_obs_var
from ._cleanup_anndata import _cleanup_anndata
from ._write_h5ad import _write_h5ad


def _compose_anndata(DataDict,
                     metadata,
                     feature_df,
                     use_matrix,
                     write_h5ad,
                     outpath,
                     silent,
                    ):
    
    adata = _dict_to_adata(DataDict)
    adata = _add_obs_var(adata, metadata, feature_df)
    adata = _cleanup_anndata(adata)
    
    if not silent:
        print(adata)
        
    if write_h5ad:
        _write_h5ad(adata, use_matrix, outpath, silent)
    
    return adata