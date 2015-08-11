from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='gsp',
      version=version,
      description="General Sequence Pattern matching",
      long_description="""\
An implementation of the GSP (general sequence pattern matching) algorithm, allowing customization and modification of most parts.""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='gsp pattern-matching',
      author='Paul Agapow',
      author_email='paul@agapow.net',
      url='http://www.agapow.net/software/gsp',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
