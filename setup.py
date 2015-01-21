#!/usr/bin/env/python

'''
Created on Jan 19, 2015

@author: sebastian
'''

from setuptools import setup

setup(
    name="rna_prediction",
    version="1.0",
    description="RNA Structure Prediction",
    author="Sebastian Ratz",
    author_email="sebastian.ratz@student.kit.edu",
    packages=["rna_prediction"],
    entry_points={
        "console_scripts": [
            "rna_predict = rna_prediction.main:main"
        ]
    },
    install_requires=["biopython >= 1.5", "numpy >= 1.6", "matplotlib > 1.1"]
)
