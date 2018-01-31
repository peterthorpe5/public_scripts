# try using distribute or setuptools or distutils.
try:
    import distribute_setup
    distribute_setup.use_setuptools()
except ImportError:
    pass

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
import sys
import re

# parse version from package/module without importing or evaluating the code
with open('pycits/__init__.py') as fh:
    for line in fh:
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group('version')
            break

if sys.version_info <= (3, 4):
    sys.stderr.write("ERROR: pycits requires Python Version 3.5 " +
                     "or above...exiting.\n")
    sys.exit(1)

setup(
    name="pycits",
    version=version,
    author="Leighton Pritchard",
    author_email="leighton.pritchard@hutton.ac.uk",
    description=''.join(["This tool enables ITS1-based classification ",
                         "for Phytophthora isolates, as part of the ",
                         "THAPBI Phyto-Threats project"]),
    license="MIT",
    keywords="genome bioinformatics sequence sequencing metabarcoding",
    platforms="Posix; MacOS X",
    url="http://widdowquinn.github.io/THAPBI-pycits",  # project home
    download_url="https://github.com/widdowquinn/THAPBI-pycits",
    scripts=['run-pycits.py'],
    packages=['pycits'],
    setup_requires=['numpy'],  # ensures numpy installation
    install_requires=['biopython',
                      'biom-format'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )
