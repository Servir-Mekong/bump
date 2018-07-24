
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='bump',
      version='0.0.1',
      description='Basic Utility Mapping Preprocessor',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/Servir-Mekong/bump',
      packages=setuptools.find_packages(),
      author='Kel Markert',
      author_email='kel.markert@nasa.gov',
      license='MIT',
      packages=['bump'],
      zip_safe=False)
