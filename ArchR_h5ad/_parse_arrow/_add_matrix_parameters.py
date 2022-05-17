
__module_name__ = "_add_matrix_parameters.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


def _return_matrix_params(file, use_matrix):

    params = file[use_matrix]["Info"]["Params"][:][0]

    params_ = []
    for val in params.tolist():
        if type(val) is bytes:
            params_.append(val.decode("utf-8"))
        else:
            params_.append(val)

    return params_

def _add_matrix_parameters(arrow, matrices=["GeneScoreMatrix", "TileMatrix"]):

    file = arrow._file

    for matrix in matrices:
        if matrix in list(file.keys()):
            key_added = "_params_{}".format(matrix)
            arrow.__setattr__(key_added, _return_matrix_params(arrow._file, matrix))