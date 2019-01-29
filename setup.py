#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(name='deepRAM',
      version='1.0.0',
      description='sequence specificities prediction of DNA- and RNA-binding proteins using deep learning approach',
      author='Ameni Trabelsi, Mohamed Chaabane and Asa Ben Hur',
      author_email='asa@cs.colostate.edu',
      maintainer='Ameni Trabelsi, Mohamed Chaabane',
      maintainer_email='chaabanemohamed2@gmail.com',
      url='https://github.com/MedChaabane/deepRAM',
      packages=find_packages(),
      install_requires=[
            'torch',
          'numpy>=1.11.2',
          'sklearn',
          'gensim',

      ],
      entry_points={
          'console_scripts': [
              'deepRAM=deepRAM:main',

          ],
      },
      )
