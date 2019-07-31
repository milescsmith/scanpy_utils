#!/usr/bin/env python
import sys
if sys.version_info < (3,):
    sys.exit('scanpy_utils requires Python >= 3.6')
from pathlib import Path
from setuptools import setup, find_packages
import versioneer

try:
    from scanpy_utils import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = ''

setup(
    name="scanpy_utils",
    version=versioneer.get_version(),
    description="utilities for working with scanpy objects",
    long_description=Path('README.rst').read_text('utf-8'),
    url="https://github.com/milescsmith/scanpy_utils",
    author=__author__,
    author_email=__email__,
    license="Proprietary",
    python_requires=">=3.6",
      install_requires=[
        l.strip() for l in
        Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    classifiers=[
        "Development Status :: 1 - Alpha",
        'Environment :: Console',
        "Intended Audience :: Developers",
        'Programming Language :: Python :: 3',
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    # extras_require=dict(
    #     doc=['sphinx', 'sphinx_rtd_theme', 'sphinx_autodoc_typehints'],
    # ),
    packages=find_packages(),
    include_package_data=True,

    keywords="scrna",
    package_dir={"scanpy_utils": "scanpy_utils"},
    # package_data={"merrycrispr": ["data/*.*"]},
    dependency_links=[
        "https://github.com/dpeerlab/Palantir/tarball/master#egg=palantir-v0.2.1",
        "https://github.com/milescsmith/scvelo/tarball/master#egg=scvelo",
        "https://github.com/jacoblevine/PhenoGraph/master#egg-phenograph",
    ],
)
