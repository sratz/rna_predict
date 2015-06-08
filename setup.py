#!/usr/bin/env/python
# coding: utf-8

from setuptools import setup


def extract_metadata():
    with open("rna_predict/__init__.py") as f:
        for line in f:
            if line.startswith('__version__'):
                exec(line.strip())
    return locals()


metadata = extract_metadata()

setup(
    name="rna_predict",
    version=metadata["__version__"],
    description="RNA Tertiary Structure Prediction",
    author="Sebastian Ratz",
    author_email="sebastian.ratz@student.kit.edu",
    packages=["rna_predict"],
    package_data={"rna_predict": ["structure_info/*"]},
    entry_points={
        "console_scripts": [
            "rna_predict = rna_predict.main:main"
        ]
    },
    install_requires=["biopython >= 1.5", "numpy >= 1.6", "matplotlib > 1.1"],
    zip_safe=False
)
