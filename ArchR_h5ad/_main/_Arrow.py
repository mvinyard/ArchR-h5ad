
__module_name__ = "_Arrow.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import h5py
import licorice_font


# import local dependencies #
# ------------------------- #
from .._parse_arrow._read_arrow_chromosome import _read_arrow_chromosome
from .._parse_arrow._add_ArchR_metadata import _add_ArchR_metadata
from .._parse_arrow._add_matrix_parameters import _add_matrix_parameters
from .._compose_adata._compose_anndata import _compose_anndata


class _Arrow:

    """Class for reading an Arrow File from .h5"""

    def __init__(
        self,
        path,
        matrices=["GeneScoreMatrix", "TileMatrix"],
        metadata_keys=["ArchRVersion", "Class"],
        silent=False,
        verbose=False
    ):

        self._path = path
        self._file = h5py.File(self._path)
        self._silent = silent
        self._verbose = verbose
        _add_ArchR_metadata(self, metadata_keys=metadata_keys)
        _add_matrix_parameters(self, matrices)

    def to_adata(self, use_matrix="GeneScoreMatrix", outpath="./", write_h5ad=True):
        
        
        self._use_matrix = use_matrix
        self._outpath = outpath
        
        if not self._silent:
            mtx = licorice_font.font_format(self._use_matrix, ["BOLD", "BLUE"])
            print("Reading ArchR {} to AnnData".format(mtx))
        
        self._DataDict = _read_arrow_chromosome(self._file, self._use_matrix, self._verbose)
        self._adata = _compose_anndata(DataDict=self._DataDict,
                                       metadata=self._file['Metadata'],
                                       feature_df=self._file[self._use_matrix]["Info"]["FeatureDF"],
                                       use_matrix=self._use_matrix,
                                       write_h5ad=write_h5ad,
                                       outpath=outpath,
                                       silent=self._silent,
                                      )
        
        