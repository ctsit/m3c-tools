###############################################################################
# Copyright 2020 University of Florida. All rights reserved.
# This file is part of UF CTS-IT's M3C project.
# Use of this source code is governed by the license found in the LICENSE file.
###############################################################################

import setuptools

VERSION = "0.13.0"

setuptools.setup(
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

    packages={"m3c": "m3c"},
    include_package_data=True,

    entry_points={
        "console_scripts": [
            "m3c = m3c.__main__:main"
        ]
    },

    install_requires=[
        "PyYAML>=5.4",
        "requests==2.23.0",
        "psycopg2-binary==2.8.5",
        "Flask==1.1.2",
        "biopython==1.76",
    ],

    python_requires=">=3.6.0",
)
