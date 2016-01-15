from distutils.core import setup
setup(
  name = 'PCRMarkerDesign',
  packages = ['PCRMarkerDesign'], # this must be the same as the name above
  version = '1.0',
  description = 'Bulk PCR Marker Design from NGS',
  author = 'John McCallum',
  author_email = 'john.mccallum@plantandfood.co.nz',
  url = 'https://github.com/cfljam/galaxy-pcr-markers', #  the URL to the github repo
  download_url = 'https://github.com/cfljam/galaxy-pcr-markers/tarball/1.0' , # Tarball URL
  keywords = ['PCR', 'genetics', 'markers','genomics'], # arbitrary keywords
  install_requires = ['cython','bcbio-gff','primer3-py'],
  scripts=['design_primers.py','vcf2gvf.py','find_CAPS.py'],
  classifiers = [],
)
