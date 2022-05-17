
__module_name__ = "_ordered_chromosomes.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import numpy as np


def _ordered_chromosomes():
    return ["chr{}".format(i) for i in np.append(np.arange(1, 23), "X")]