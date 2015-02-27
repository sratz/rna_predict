#!/usr/bin/env/python
# coding: utf-8

from setuptools import setup


def extract_metadata():
    with open("rna_prediction/__init__.py") as f:
        for line in f:
            if line.startswith('__version__'):
                exec(line.strip())
    return locals()


metadata = extract_metadata()

setup(
    name="rna_prediction",
    version=metadata["__version__"],
    description="RNA Tertiary Structure Prediction",
    author="Sebastian Ratz",
    author_email="sebastian.ratz@student.kit.edu",
    packages=["rna_prediction"],
    package_data={"rna_prediction": ["structure_info/*"]},
    entry_points={
        "console_scripts": [
            "rna_predict = rna_prediction.main:main"
        ]
    },
    install_requires=["biopython >= 1.5", "numpy >= 1.6", "matplotlib > 1.1"],
    zip_safe=False
)
