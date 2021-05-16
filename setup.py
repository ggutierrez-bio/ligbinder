#!/usr/bin/env python

from setuptools import setup

setup(
    name="LigBinder",
    version="0.1.4",
    description="automatic targeted molecular dynamics for ligand bindning",
    author="Guillermo Gutierrez",
    author_email="",
    url="https://github.com/ggutierrez-bio/ligbinder",
    packages=["ligbinder"],
    inlcude_package_data=True,
    data_files=[("ligbinder", ["ligbinder/data/default_config.yml"])],
    scripts=["bin/ligbinder"],
    install_requires=["numpy", "pyyaml", "pytraj>=2", "parmed"]
)
