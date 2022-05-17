
__module_name__ = "_read_ArchR_to_adata.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import local dependencies #
# ------------------------- #
from ._Arrow import _Arrow


def _read_arrow_to_adata(
    path,
    matrices=["GeneScoreMatrix", "TileMatrix"],
    metadata_keys=["ArchRVersion", "Class"],
    use_matrix="GeneScoreMatrix",
    silent=False,
    write_h5ad=True,
):
    """
    Read an ArchR ".arrow" file as AnnData (adata).
 
    Parameters:
    -----------
    path
        path to ArchR.arrow file.
        type: str

    matrices
        Matrices saved in the ArchR.arrow file.
        default: ["GeneScoreMatrix", "TileMatrix"]
        type: list(str)

    metadata_keys
        Keys to high-level metadata saved by ArchR.
        default: ["ArchRVersion", "Class"]
        type: list(str)

    use_matrix
        Which matrix to use. Currently TileMatrix is not implemented.
        default: "GeneScoreMatrix",
        type: str

    silent
        If True, print extra messages.
        default: False
        type: bool
        
    write_h5ad
        default: True
        type: bool
    
    Returns:
    --------
    adata
        anndata._core.anndata.AnnData

    Notes:
    ------
    (1)
    """

    arrow = _Arrow(
        path, matrices=matrices, metadata_keys=metadata_keys, silent=silent
    )
    arrow.to_adata(use_matrix=use_matrix, write_h5ad=write_h5ad)

    return arrow._adata