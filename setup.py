from distutils.core import setup
from ductape import __version__ 

setup(
    name = 'DuctApe',
    version = __version__,
    author = 'Marco Galardini',
    author_email = 'marco.galardini@unifi.it',
    packages = ['ductape','ductape.common', 'ductape.genome', 'ductape.kegg',
                'ductape.phenome', 'ductape.storage', 'ductape.storage.SQLite',
				'ductape.storage.data'],
	package_data={'ductape.storage.data': ['*.tsv']},
    scripts = ['dape', 'dgenome', 'dphenome'],  
    url = 'https://combogenomics.github.com/DuctApe', 
    license = 'LICENSE.txt',
    description = 'Analyzing and linking genomics and phenomics experiments',
    long_description = open('README.txt').read(),
    install_requires = ['argparse >= 1.1', 'biopython >= 1.5', 'numpy',
                        'matplotlib >= 1.1', 'multiprocessing',
                        'scipy', 'scikit-learn >= 0.11', 'PyYAML',
			'networkx']
)



