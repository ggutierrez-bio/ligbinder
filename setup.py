#!/usr/bin/env python

from setuptools import setup

setup(
    name='LigBinder',
    version='0.1',
    description='automatic targeted molecular dynamics for ligand bindning',
    author='Guillermo Gutierrez',
    author_email='',
    url='',
    packages=['ligbinder'],
    # as an example of additional data:
    # inlcude_package_data=True,
    # package_data={'ligbinder.data':[ligbinder/data/*]}
    scripts=['ligbinder'],
    install_requires=["pytraj>=2"]
)
