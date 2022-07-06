from setuptools import setup
import re
import os
import sys


setup(
    name="ArchR_h5ad",
    version="0.0.12",
    python_requires=">3.6.0",
    author="Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard",
    author_email="mvinyard@broadinstitute.org",
    url="",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="ArchR_h5ad: Read .arrow files (from ArchR) to anndata.",
    packages=[
        "ArchR_h5ad",
        "ArchR_h5ad._compose_adata",
        "ArchR_h5ad._main",
        "ArchR_h5ad._parse_arrow",
        "ArchR_h5ad._utility_functions",
    ],
    install_requires=[
        "anndata>=0.7.8",
        "licorice_font>=0.0.3",
        "tqdm>=4.64.0",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
