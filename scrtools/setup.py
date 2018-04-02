from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()

requires = ['scanpy>=0.4',
			'anndata>=0.5',
			'matplotlib>=2.0.0',
			'pandas>=0.21',
			'scipy',
			'seaborn',
			'psutil',
			'h5py',
			'xlrd',
			'scikit-learn>=0.19.1',
			'statsmodels',
			'networkx',
			'natsort',
			'joblib',
			'profilehooks',
			'numpy',
			'tables',
			'xlsxwriter',
			'fisher',
			'tqdm',
			'loompy',
			'louvain',
			'MulticoreTSNE',
			'docopt',
			'setuptools'
		  ]

setup(
	name='scrtools',
	version='0.1.0',
	description='scRNA-Seq analysis tools that scale to millions of cells, built upon scanpy',
	long_description = long_description,
	url='https://github.com/broadinstitute/scRNA-Seq/tree/master/scrtools',
	author='Bo Li, Joshua Gould',
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
		'scrtools.annotate_cluster': ['cell_type_markers.json'],
	},
	entry_points={
		'console_scripts': [
			'scrtools=scrtools.__main__:main',
		],
	},
)
