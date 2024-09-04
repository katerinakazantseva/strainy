import os
import sys
import subprocess
import shutil

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

from strainy.__version__ import __version__


setup(name='strainy',
      version=__version__,
      description='Metagenomic strain phasing and assmebly using long reads',
      url='https://github.com/katerinakazantseva/strainy',
      author='Ekaterina Kazantseva & Ataberk Donmez',
      author_email = 'ekaterina.v.kazantseva@gmail.com',
      license='CC BY-NC-SA 4.0',
      packages=['strainy', 'strainy/clustering', 'strainy/gfa_operations', 'strainy/reports', 'strainy/simplification'],
      entry_points={'console_scripts': ['strainy = strainy.main:main']},
      )
