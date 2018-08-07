# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

exec(open('nanoplotter/version.py').read())
here = path.abspath(path.dirname(__file__))


setup(
    name='nanoplotter',
    version=__version__,
    description='Plotting functions of Oxford Nanopore sequencing data',
    long_description=open(path.join(here, "README.md")).read(),
    long_description_content_type="text/markdown",
    url='https://github.com/wdecoster/nanoplotter',
    author='Wouter De Coster',
    author_email='decosterwouter@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='nanopore sequencing plotting quality control',
    packages=find_packages() + ['scripts'],
    python_requires='>=3',
    setup_requires=['setuptools>=38.6.0'],
    install_requires=['pandas',
                      'numpy',
                      'scipy',
                      'matplotlib>=2.1.0',
                      'seaborn>=0.9.0',
                      "pauvre==0.1.86",
                      "statsmodels>=0.8.0",
                      'plotly>=2.5.1'],
    package_dir={'nanoplotter': 'nanoplotter'})
