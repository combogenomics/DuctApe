from distutils.core import setup
from ductape import __version__ 
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md')) as f:
        long_description = f.read()

setup(
    name = 'DuctApe',
    version = __version__,
    author = 'Marco Galardini',
    author_email = 'mgala@bu.edu',
    packages = ['ductape','ductape.common', 'ductape.genome', 'ductape.kegg',
                'ductape.phenome', 'ductape.storage', 'ductape.storage.SQLite',
				'ductape.storage.data'],
	package_data={'ductape.storage.data': ['*.tsv']},
    scripts = ['dape', 'dgenome', 'dphenome'],  
    url = 'https://combogenomics.github.com/DuctApe', 
    license = 'LICENSE.txt',
    description = 'Analyzing and linking genomics and phenomics experiments',
    long_description = long_description,
    long_description_content_type='text/markdown',
    install_requires = ['argparse >= 1.1', 'biopython >= 1.5', 'numpy',
                        'matplotlib >= 1.1',
                        'scipy', 'scikit-learn >= 0.11', 'PyYAML',
			'networkx']
)



