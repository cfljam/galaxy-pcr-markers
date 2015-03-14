[![Stories in Ready](https://badge.waffle.io/cfljam/galaxy-pcr-markers.png?label=ready&title=Ready)](https://waffle.io/cfljam/galaxy-pcr-markers)
galaxy-pcr-markers
==================

Scripts for design of PCR-based Marker Assays from DNA sequence variant data. 

Xml wrappers are included for use in the  Galaxy  workflow environment.
Also available for download at Galaxy Toolshed http://toolshed.g2.bx.psu.edu/
hg clone http://toolshed.g2.bx.psu.edu/repos/john-mccallum/pcr_markers. 

NOTE that the primer design tool *design_primers.py*  no longer relies on EMBOSS ePrimer3 (which depends on Primer3 1.1.4.)
 and should work with Primer3 V2.3+

This release adds optimised design of high-resolution melting PCR assays using the uMelt web service provided by the Wittwer Lab at University of Utah https://www.dna.utah.edu/umelt/umelt.html


**CITATION**
A Toolkit For Bulk PCR-Based Marker Design From Next-Generation Sequence Data: Application For Development Of A Framework Linkage Map In Bulb Onion (Allium cepa L.) (2012)

Samantha Baldwin, Roopashree Revanna, Susan Thomson, Meeghan Pither-Joyce, Kathryn Wright, Ross Crowhurst, Mark Fiers, Leshi Chen, Richard MacKnight, John A. McCallum

BMC Genomics 2012, 13:637  http://www.biomedcentral.com/1471-2164/13/637/abstract

uMELT: prediction of high-resolution melting curves and dynamic melting profiles of PCR products in a rich web application.
Zachary Dwight1, Robert Palais and Carl T. Wittwer http://bioinformatics.oxfordjournals.org/content/27/7/1019

**Acknowledgements**
Development of these tools was funded by the New Zealand Ministry for Business, Innovation & Employment project 'Virtual Institute of Statistical Genetics' (VISG)
See http://www.visg.co.nz
