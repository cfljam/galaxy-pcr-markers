[![Stories in Ready](https://badge.waffle.io/cfljam/galaxy-pcr-markers.png?label=ready&title=Ready)](https://waffle.io/cfljam/galaxy-pcr-markers)
galaxy-pcr-markers
==================

Scripts for design of PCR-based Marker Assays from DNA sequence variant data and optimised design of high-resolution melting PCR assays using the uMelt web service provided by the Wittwer Lab at University of Utah https://www.dna.utah.edu/umelt/umelt.html

Xml wrappers are included for use in the  Galaxy  workflow environment.
Also available for download at Galaxy Toolshed http://toolshed.g2.bx.psu.edu/
hg clone http://toolshed.g2.bx.psu.edu/repos/john-mccallum/pcr_markers.

The primer design tool *design_primers.py*  now uses the excellent [primer3-py](https://github.com/benpruitt/primer3-py) See  http://benpruitt.github.io/primer3-py/index.html

Dependencies
------------
- Linux or OSX - not tested on Windows 
- Python 2.7
- NumPy
- BioPython see https://pypi.python.org/pypi/biopython
- [BcBio-gff](https://github.com/chapmanb/bcbb/tree/master/gff)
- [Primer3-py](https://github.com/benpruitt/primer3-py)

We recommend running inside a Conda environment

- install [Miniconda](https://conda.io/docs/install/quick.html)
- create a fresh environment and activate it

```
conda create -y -n Py2PCR python=2.7
source activate Py2PCR
```
Install the dependencies using pip

```
pip install numpy
pip install biopython
pip install bcbio-gff primer3-py
```
N.B. Primer3 install is **NOT** required now since design is handled by *primer3-py*


**CITATION**
A Toolkit For Bulk PCR-Based Marker Design From Next-Generation Sequence Data: Application For Development Of A Framework Linkage Map In Bulb Onion (Allium cepa L.) (2012)

Samantha Baldwin, Roopashree Revanna, Susan Thomson, Meeghan Pither-Joyce, Kathryn Wright, Ross Crowhurst, Mark Fiers, Leshi Chen, Richard MacKnight, John A. McCallum

BMC Genomics 2012, 13:637  http://www.biomedcentral.com/1471-2164/13/637/abstract

uMELT: prediction of high-resolution melting curves and dynamic melting profiles of PCR products in a rich web application.
Zachary Dwight1, Robert Palais and Carl T. Wittwer http://bioinformatics.oxfordjournals.org/content/27/7/1019

**Acknowledgements**
Development of these tools was funded by the New Zealand Ministry for Business, Innovation & Employment project 'Virtual Institute of Statistical Genetics' (VISG)
See http://www.visg.co.nz
