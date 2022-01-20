from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from codecs import open
from pathlib import Path
from os import path


here = path.abspath(path.dirname(__file__))
with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

extensions = [
    Extension("pegasus.cylib.fast_utils", ["ext_modules/fast_utils.pyx"]),
    Extension("pegasus.cylib.cfisher", ["ext_modules/cfisher.pyx"]),
    Extension("pegasus.cylib.de_utils", ["ext_modules/diff_expr_utils.pyx"]),
]

setup(
    name="pegasuspy",
    use_scm_version=True,
    description="Pegasus is a Python package for analyzing sc/snRNA-seq data of millions of cells",
    long_description=long_description,
    url="https://github.com/klarman-cell-observatory/pegasus",
    author="Yiming Yang, Joshua Gould and Bo Li",
    author_email="cumulus-support@googlegroups.com",
    classifiers=[  # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "Framework :: Jupyter",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Topic :: Software Development :: Build Tools",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords="single cell/nucleus genomics analysis",
    packages=find_packages(),
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    ext_modules=cythonize(extensions),
    setup_requires=["Cython", "setuptools_scm"],
    extras_require=dict(
        tsne=["fitsne"],
        louvain=["louvain"],
        scanorama=["scanorama"],
        torch=["torch", "harmony-pytorch", "nmf-torch"],
        forceatlas=["forceatlas2-python"],
        mkl=["mkl"],
        all=["fitsne", "louvain", "scanorama", "torch", "harmony-pytorch", "nmf-torch", "forceatlas2-python"]
    ),
    python_requires="~=3.6",
    package_data={
        "pegasus.annotate_cluster": [
            "human_immune_cell_markers.json",
            "mouse_immune_cell_markers.json",
            "mouse_brain_cell_markers.json",
            "human_brain_cell_markers.json",
            "human_lung_cell_markers.json",
        ],
        "pegasus.check_sample_indexes": ["chromium-shared-sample-indexes-plate.json", "Chromium-i7-Multiplex-Kit-N-Set-A-sample-indexes-plate.json"],
        "pegasus": ["data_files/*.gmt"],
    },
    entry_points={"console_scripts": ["pegasus=pegasus.__main__:main"]},
)
