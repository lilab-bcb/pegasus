from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()

requires = ['anndata==0.6.4',
			'matplotlib>=2.0.0',
			'pandas>=0.21',
			'Cython',
			'scipy',
			'seaborn',
			'scikit-learn>=0.19.1',
			'statsmodels',
			'natsort',
			'numpy',
			'tables',
			'xlsxwriter',
			'fisher',
			'loompy',
			'louvain',
			'MulticoreTSNE',
			'docopt',
			'setuptools',
			'plotly',
			'pybind11',
			'umap-learn',
			'fitsne',
			'hdbscan',
			'pyarrow',
			'google-cloud-bigquery'
			# 'hnswlib'
		  ]

setup(
	name='scCloud',
	version='0.6.0',
	description='scRNA-Seq analysis tools that scale to millions of cells',
	long_description = long_description,
	url='https://github.com/broadinstitute/scRNA-Seq/tree/master/scCloud',
	author='Bo Li, Joshua Gould, Siranush Sarkizova, Marcin Tabaka, Orr Ashenberg, and et al.',
	author_email='libo@broadinstitute.org',
	classifiers=[ # https://pypi.python.org/pypi?%3Aaction=list_classifiers
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Developers',
		'Intended Audience :: Science/Research',
		'Topic :: Software Development :: Build Tools',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Programming Language :: Python :: 3',
		'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 3.6',
		'Framework :: Jupyter',
		'Natural Language :: English',
		'Operating System :: MacOS :: MacOS X',
		'Operating System :: POSIX :: Linux',
		'Topic :: Scientific/Engineering :: Bio-Informatics'
	],
	keywords='single cell RNA-Seq data analysis', 
	packages=find_packages(),
	install_requires=requires,
	python_requires='~=3.5',
	package_data={
		'scCloud.annotate_cluster': ['human_immune_cell_markers.json', 'mouse_immune_cell_markers.json', 'mouse_brain_cell_markers.json', 'human_brain_cell_markers.json'],
		'scCloud': ['ext/GraphLayout.java', 'ext/GraphLayout.class', 'ext/gephi-toolkit-0.9.2-all.jar']
	},
	entry_points={
		'console_scripts': [
			'scCloud=scCloud.__main__:main',
		],
	},
)
