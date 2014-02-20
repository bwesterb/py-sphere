#!/usr/bin/env python

import sys
from setuptools import setup
from src.version import get_version

install_requires = [
        ]

setup(
    name='sphere',
    version=get_version(),
    description='Geometry on the sphere',
    author='Bas Westerbaan',
    author_email='bas@westerbaan.name',
    url='https://github.com/bwesterb/py-sphere',
    packages=['sphere',
              'sphere.version',
              'sphere.tests'],
    package_dir={'sphere': 'src'},
    license='GPL 3.0',
    install_requires=install_requires,
    # TODO add classifiers
    classifiers = [
        ],
    test_suite='sphere.tests',
    ),
