#!/usr/bin/env python3
import sys
from distutils.core import setup
from setuptools import setup, find_packages

CURRENT_PYTHON = sys.version_info[:2]
REQUIRED_PYTHON = (3, 3)

# This check and everything above must remain compatible with Python 2.7.
if CURRENT_PYTHON < REQUIRED_PYTHON:
    sys.stderr.write("graphdraw requires Python 3.3 or higher and "
                     "you current version is {}".format(CURRENT_PYTHON))
    sys.exit(1)

reqs = []
with open("requirements.txt", "r") as f:
    for l in f:
        reqs.append(l.strip())

setup(name='graphdraw',
      version='0.0.1',
      description='For manipulation and drawing of graphs',
      author='Fawaz Dabbaghie',
      author_email='fawaz@hhu.de',
      url='https://github.com/marschall-lab/project-ccl-seq',
      packages=find_packages(),
      # scripts=['bin/main.py'],
      license="LICENSE",
      long_description=open("README.md").read(),
      long_description_content_type='text/markdown',
      install_requires=reqs,
      include_package_data=True,
      entry_points={
          "console_scripts": [
              "graphdraw=graphdraw.main:main"
          ]}
      )
