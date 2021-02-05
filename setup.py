#!/usr/bin/env python

from setuptools import setup

setup(
    name='Samplematic - BigChem',
    version='0.1',
    description='automatic targeted molecular dynamics',
    author='Guillermo Gutierrez',
    author_email='',
    url='',
    packages=['samplematic-bigchem'],
    # as an example of additional data:
    # inlcude_package_data=True,
    # package_data={'samplematic.data':[samplematic/data/*]}
    scripts=['samplematic'],
    install_requires=["pytraj>=2"]
)
