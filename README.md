DuctApe
=======

Analyzing and linking genomics and phenomics experiments
--------------------------------------------------------

DuctApe is a python library with a series of command-line programs that will help bioinformaticians to analyze genomes AND phenomic experiments.
The final purpose of the program is to combine the genomic informations (encoded as KEGG pathways) with the results of phenomic experiments (for now the Phenotype microarrays are supported) and highlight the genes that may be responsible for phenotypic variations.
Several scientific python libraries are used to perform this tasks: all the data is stored in a portable SQLite database to allow an easy exchange of results.
The most complex computing tasks are speeded up by parallelization.
There are three distinct programs with git-like syntax to ensure flexibility in the analysis.

The program is currently in the alpha phase, still lacking the utilities to combine the genomic and phenomic analyses.
Future releases will also include GUIs and web interfaces.

Requirements
------------
* Biopython
* Numpy
* SciPy
* matplotlib
* scikits.learn
* SOAPpy
* multiprocessing
* (networkx)

Installation
------------
* sudo python setup.py install

* python setup.py sdist
* sudo pip install dist/DuctApe-X.X.X.tar.gz

Quick how-to
------------
Three command line utilities will be installed

* dape
    * Used to initialize the project and add/remove the organism
    * (Performs the combination of genomics and phenomics experiments)

* dgenome
    * Handles protein sequences for each organism of the project
    * For pangenomic studies it can build the pangenome using the BBH algorithm (user-defined pangenomes can be added)
    * The proteins are mapped to KEGG metabolic pathways

* dphenome
    * Handles phenomic experiments for each organism of the project (Phenotype microarray experiments)
    * Performs zero-subtraction on the biolog plates
    * Calculates the growth parameters
    * Ranks the single experiments with a clusterization over growth parameters
    * The compounds from the phenomic experiments are mapped to KEGG

Examples
--------
* Single organism experiment
    * dape init (initializes the project)
    * dape add MyOrg (adds my organism using the ID MyOrg)

    * dgenome add MyOrg.faa MyOrg (adds the proteome of MyOrg)
    * dgenome add-ko MyOrg.tab (adds the output of KAAS, a KEGG mapper)
    * dgenome start (maps the proteome to KEGG)
    * dgenome map (outputs the KEGG metabolic maps)
    * dgenome stats (statistic and graphics)
    * dgenome export (exports the genomic data)

    * dphenome add MyOrg.csv MyOrg (adds the phenomic experiment, BIOLOG data)
    * dphenome zero (performs control subtraction)
    * dphenome start -n 4 (calculates the growth parameters and performs the clusterization, using 4 CPUs)
    * dphenome plot (plots the growth curves)
    * dphenome purge -d 3 keep-max (removes inconsistent replicas: keep the highest replicas when there is an activity index delta >= 3)
    * dphenome plot (plots only those curves that are not purged)
    * dphenome restore (restore the purged replicas)
    * dphenome stats
    * dphenome export

    * (dape map)

* Mutant experiment
    * dape init
    * dape add MyOrg
    * dape add-mut -m MyOrg -k deletion MyMut (adds mutant MyMut, a deletion mutant of MyOrg)

    * dgenome add-dir MyFolder (adds the proteins files found in this directory)
    * dgenome add-ko MyOrg.tab
    * dgenome start
    * dgenome map MyMut (plots only tha maps of the mutant)

    * dphenome add-dir MyPhenomicFolder (adds the phenomic files found in this directory)
    * dphenome zero
    * dphenome start -n 4
    * dphenome purge -d 3 keep-max
    * dphenome plot

    * (dape map)

* Pangenomic experiment
    * dape init
    * dape add MyOrg
    * dape add MyOrg2
    * dape add MyOrg3

    * dgenome add-dir MyFolder
    * dgenome add-ko MyOrg.tab MyOrg2.tab MyOrg3.tab
    * dgenome start -n 4 (also performs pangenome creation using 4 CPUs)
    * dgenome map (plots the maps for the whole pangenome)

    * dphenome add-dir MyPhenomicFolder
    * dphenome zero
    * dphenome start -n 4
    * dphenome purge -d 3 keep-max
    * dphenome plot

    * (dape map)

More informations
-----------------
Each program options and parameters can be queried adding -h

Contacts
--------
This program has been developed by Marco Galardini, Department of evolutionary biology, University of Florence
