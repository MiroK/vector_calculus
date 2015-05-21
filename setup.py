#!/usr/bin/env python

from distutils.core import setup
setup(name = 'vector_calculus',
      version = '0.0.1',
      description = 'Very restricted vector calculus in cartesian basis',
      author = 'Miroslav Kuchta',
      author_email = 'mirok@math.uio.no',
      url = '',
      packages = ['vector_calculus',
                  'vector_calculus.containers',
                  'vector_calculus.operators']
)
