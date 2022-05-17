
__module_name__ = "_h5ad_filename.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import licorice_font
import os


def _h5ad_filepath(adata, use_matrix, outpath):
    
    filename = "{}.{}.h5ad".format(adata.uns['metadata_dict']['Sample'], use_matrix)
    return os.path.join(outpath, filename)

def _write_h5ad(adata, use_matrix, outpath, silent):

    h5ad_filepath = _h5ad_filepath(adata, use_matrix, outpath)
    if not silent:
        msg = licorice_font.font_format("Saving to", ["BOLD"])
        print("\n{}: {}".format(msg, h5ad_filepath))

    adata.write_h5ad(h5ad_filepath)