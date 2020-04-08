###############################################################################
# Copyright 2020 University of Florida. All rights reserved.
# This file is part of UF CTS-IT's M3C project.
# Use of this source code is governed by the license found in the LICENSE file.
###############################################################################

from setuptools import setup, find_packages

VERSION="0.7.1"

setup(
    name="m3c",
    version=VERSION,
    author="UF CTS-IT",
    author_email="ctsit@ctsi.ufl.edu",
    maintainer="UF CTS-IT",
    maintainer_email="ctsit@ctsi.ufl.edu",
    url="https://github.com/ctsit/metab_import",
    license="Apache 2.0",
    description="Collection of M3C tools",
    keywords=["metabolomics", "tpf", "vivo"],
    download_url="https://github.com/ctsit/m3c/releases/tag/" + VERSION,

    package_dir = {'m3c': 'm3c'},
    packages = find_packages(),

    entry_points={
        "console_scripts": [
            "m3c = m3c.__main__:main"
        ]
    },

    install_requires=[
        "PyYAML==5.1.2",
        "requests==2.22.0",
        "psycopg2-binary==2.8.3",
        "flask==1.1.1",
        "biopython==1.74",
    ],

    python_requires=">=3.6.0",
)
