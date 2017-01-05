__author__ = 'Sergei F. Kliver'

from setuptools import setup, find_packages
from os.path import join, dirname


setup(name='KrATER',
      version='0.17',
      packages=find_packages(),
      author='Sergei F. Kliver',
      url='https://github.com/mahajrod/KRATER',
      author_email='mahajrod@gmail.com',
      install_requires=['scipy', 'numpy', 'matplotlib'],
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      scripts=['draw_kmer_distribution_from_fastq.py',
               'draw_kmer_distribution_from_histo.py',
               'draw_kmer_distribution_from_jellyfish_database.py'])
