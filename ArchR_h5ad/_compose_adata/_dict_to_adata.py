
__module_name__ = "_dict_to_adata.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import anndata
import scipy.sparse


# import local dependencies #
# --------------- #
from .._utility_functions._ordered_chromosomes import _ordered_chromosomes


def _dict_to_adata(DataDict):

    _ordered_matrices = []
    for chrom in _ordered_chromosomes():
        if not DataDict[chrom] is None:
            _ordered_matrices.append(DataDict[chrom])
            
    X_ = scipy.sparse.hstack(_ordered_matrices)
    
    return anndata.AnnData(X_, dtype=X_.dtype)