[![Stories in Ready](https://badge.waffle.io/cfljam/galaxy-pcr-markers.png?label=ready&title=Ready)](https://waffle.io/cfljam/galaxy-pcr-markers)
galaxy-pcr-markers
==================

Scripts for design of PCR-based Marker Assays from DNA sequence variant data and optimised design of high-resolution melting PCR assays using the uMelt web service provided by the Wittwer Lab at University of Utah https://www.dna.utah.edu/umelt/umelt.html

Xml wrappers for use in the  Galaxy  workflow environment are deprecated and not maintained. 
See older versions 

We are currently (2017)  refactoring to enable use with VCF and BED formats on genome-scale data.  

The primer design tool *design_primers.py*  now uses the excellent [primer3-py](https://github.com/benpruitt/primer3-py) See  http://benpruitt.github.io/primer3-py/index.html

Dependencies
------------
- Linux or OSX - not tested on Windows 
- Python 2.7
- NumPy
- SciPy
- BioPython see https://pypi.python.org/pypi/biopython
- [BcBio-gff](https://github.com/chapmanb/bcbb/tree/master/gff)
- [Primer3-py](https://github.com/benpruitt/primer3-py)

Recommended
-----------

- ipython or jupyter

Installation
-----------

We recommend running inside a Conda environment

- install [Miniconda](https://conda.io/docs/install/quick.html)
- create a fresh environment and activate it

```
conda create -y -n Py2PCR python=2.7
source activate Py2PCR
```
### Install the dependencies using pip

```
pip install numpy scipy
pip install biopython
pip install bcbio-gff primer3-py
pip install ipython==4.2
```
N.B. Primer3 install is **NOT** required now since design is handled by *primer3-py*

### Clone or download the repo and move into it

```
git clone https://github.com/cfljam/galaxy-pcr-markers
cd galaxy-pcr-markers
```
### Check all is well by running on small test data, specifying one primer set and melt prediction for HRM

>python design_primers.py -i ./test-data/targets.fasta -g  ./test-data/targets.gff -T ./test-data/targets -n 1 -u

It should return:
```
SNP_Target_ID Position Ref_base Variant_base Amplicon_bp PRIMER_LEFT_SEQUENCE PRIMER_RIGHT_SEQUENCE ref_melt_Tm var_melt_Tm Tm_difference
k69_93535:SAMTOOLS:SNP:1147 1147 C G 285 CTCTTCAGTTGCTTCCTGCC CTTCACTCCTTCTCGCGTTC 87.35 87.7 0.35
k69_93535:SAMTOOLS:SNP:1336 1336 G A 149 GAACGCGAGAAGGAGTGAAG GCAACCCAGGTTTCAACTCC 88.75 88.7 0.05
k69_98089:SAMTOOLS:SNP:550 550 G A 227 GGAGAAGGTCGAGGTCAGC ACGGCCGAATATACATACAACG 85.75 86.2 0.45
k69_98089:SAMTOOLS:SNP:625 625 A G 227 GGAGAAGGTCGAGGTCAGC ACGGCCGAATATACATACAACG 85.75 86.25 0.5

------------------------------

**CITATION**
A Toolkit For Bulk PCR-Based Marker Design From Next-Generation Sequence Data: Application For Development Of A Framework Linkage Map In Bulb Onion (Allium cepa L.) (2012)

Samantha Baldwin, Roopashree Revanna, Susan Thomson, Meeghan Pither-Joyce, Kathryn Wright, Ross Crowhurst, Mark Fiers, Leshi Chen, Richard MacKnight, John A. McCallum

BMC Genomics 2012, 13:637  http://www.biomedcentral.com/1471-2164/13/637/abstract

uMELT: prediction of high-resolution melting curves and dynamic melting profiles of PCR products in a rich web application.
Zachary Dwight1, Robert Palais and Carl T. Wittwer http://bioinformatics.oxfordjournals.org/content/27/7/1019

**Acknowledgements**
Development of these tools was funded by the New Zealand Ministry for Business, Innovation & Employment project 'Virtual Institute of Statistical Genetics' (VISG)
See http://www.visg.co.nz
