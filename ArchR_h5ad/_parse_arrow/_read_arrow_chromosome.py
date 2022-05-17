
__module_name__ = "_read_arrow_chromosome.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import numpy as np
import scipy.sparse
from tqdm.notebook import tqdm


# import local dependencies #
# ------------------------- #
from .._utility_functions._ordered_chromosomes import _ordered_chromosomes


def _return_sum_chromosome_axis_sums(chromosome):

    colsums = np.array(chromosome["colSums"]).flatten().sum()
    rowsums = np.array(chromosome["rowSums"]).flatten().sum()

    return np.array([colsums, rowsums]).sum()

def _return_zero_rows(colsums):
    return np.where(colsums == 0)[0]

def _get_matrix_size(chromosome):

    """"""

    ncols = np.array(chromosome["colSums"]).flatten().shape[0]
    nrows = np.array(chromosome["rowSums"]).flatten().shape[0]

    return [ncols, nrows]

def _initialize_empty_chromosome_data_matrix(chromosome):

    [ncols, nrows] = _get_matrix_size(chromosome)

    return np.zeros([ncols, nrows])

def _return_jLengths(chromosome):
    return np.append(0, np.array(chromosome["jLengths"]).flatten()).cumsum()

def _fetch_chromosome_data_from_arrow_h5(chromosome, binary):

    colsums = np.array(chromosome["colSums"]).flatten()
    rowsums = np.array(chromosome["rowSums"]).flatten()

    axis_sums = _return_sum_chromosome_axis_sums(chromosome)
    zero_rows = _return_zero_rows(colsums)

    if axis_sums.sum() == 0:
        print("\tNo features / cells found in chromosome...")
        return None

    else:

        X_empty = _initialize_empty_chromosome_data_matrix(chromosome)
        j_lengths = _return_jLengths(chromosome)
        
        i = np.array(chromosome["i"]).flatten()
        
        if not binary:
            x = np.array(chromosome["x"]).flatten()
        else:
            x = np.ones(len(i))

        row_adj = 0
        row_sums = []

        for row in range(len(X_empty)):
            if not row in zero_rows:
                j_len_i = j_lengths[row_adj]
                if not row_adj == len(j_lengths):
                    j_len_j = j_lengths[int(row_adj + 1)]
                else:
                    j_len_j = j_len_i
                row_vals = x[j_len_i:j_len_j]
                row_sums.append(row_vals.sum())
                idx = i[j_len_i:j_len_j] - 1
                row_adj += 1
                X_empty[row, idx] = row_vals

        return scipy.sparse.csr_matrix(X_empty)
    

def _read_arrow_chromosome(h5_file, use_matrix="GeneScoreMatrix", verbose=False):

    chromosomes = list(h5_file[use_matrix].keys())
    chromosomes.remove("Info")
    
    if use_matrix == "TileMatrix":
        binary = True
    else:
        binary = False

    DataDict = {}
    if verbose:
        print("Loading chromosomes from Arrow:")
    for chrom_key in tqdm(_ordered_chromosomes(), desc="Chromosomes"):
        if chrom_key in chromosomes:
            chromosome = h5_file[use_matrix][chrom_key]
            if verbose:
                print("- {}".format(chrom_key))
            DataDict[chrom_key] = _fetch_chromosome_data_from_arrow_h5(chromosome, binary)
        else:
            print(" - Warning: {} not detected!".format(chrom_key))
    return DataDict