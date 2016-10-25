__author__ = 'mahajrod'

from setuptools import setup, find_packages
from os.path import join, dirname


setup(name='KrATER',
      version='0.1',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=['scipy', 'numpy', 'matplotlib'],
      long_description=open(join(dirname(__file__), 'README.md')).read(),)
