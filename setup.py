import versioneer
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

requires = [
    "matplotlib>=2.0.0",
    "pandas>=0.21",
    "Cython",
    "scipy<1.3",
    "seaborn",
    "scikit-learn>=0.21.3",
    "statsmodels",
    "natsort",
    "anndata",
    "numba",
    "numpy",
    "tables",
    "xlsxwriter",
    "loompy",
    "docopt",
    "setuptools",
    "plotly",
    "pybind11",
    "joblib",
    "scikit-misc",
    "pyarrow",
    "umap-learn>=0.3.9",
    "lightgbm==2.2.1",
    "python-igraph",
    "MulticoreTSNE-modified",
    "hnswlib",
    "fisher-modified",
    "louvain-github",
    "leidenalg",
    "forceatlas2-python",
    "scplot"
]

setup(
    name="sccloud",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="scRNA-Seq analysis tools that scale to millions of cells",
    long_description=long_description,
    url="https://github.com/klarman-cell-observatory/scCloudPy",
    author="Bo Li, Joshua Gould, Yiming Yang, Siranush Sarkizova",
    author_email="sccloud@googlegroups.com, sccloud@broadinstitute.org",
    classifiers=[  # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Framework :: Jupyter",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="single cell/nucleus genomics analysis",
    packages=find_packages(),
    install_requires=requires,
    extras_require=dict(fitsne=["fitsne"], mkl=["mkl"]),
    python_requires="~=3.5",
    package_data={
        "sccloud.annotate_cluster": [
            "human_immune_cell_markers.json",
            "mouse_immune_cell_markers.json",
            "mouse_brain_cell_markers.json",
            "human_brain_cell_markers.json",
        ],
        "scCloud.check_sample_indexes": ["chromium-dna-sample-indexes-plate.json"],
    },
    entry_points={"console_scripts": ["sccloud=sccloud.__main__:main"]},
)
