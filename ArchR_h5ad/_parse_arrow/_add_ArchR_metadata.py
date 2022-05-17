
__module_name__ = "_add_ArchR_metadata.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


def _add_ArchR_metadata(arrow, metadata_keys=["ArchRVersion", "Class"]):
    for key in metadata_keys:
        arrow.__setattr__("_{}".format(key), arrow._file[key][0].decode("utf-8"))