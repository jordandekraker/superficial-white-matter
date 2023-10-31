#!/usr/bin/env python

"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
"""

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="swm",
    version="0.1.0",
    author="Jordan DeKraker",
    author_email="jordandekraker@gmail.com",
    description="translates white matter surfaces inward along a Laplace gradient",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/jordandekraker/superficial-white-matter",
    packages=setuptools.find_packages(),
    license="BSD 3-Clause License",
    package_data={
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords='hippunfold registration',
    python_requires=">=3.7.0",
    install_requires=[
        "Path>=16.4.0",
        "nibabel>=3.2.2",
        "numpy>=1.16.5",
        "pandas>=0.23",
        "scipy>=1.3.3",
        "astropy>=4.3.1",
        "scikit-fmm>=2022.04.02",
    ],
    extras_require={"dev": ["gitpython", "hcp-utils", "mypy", "plotly", "pytest"]},
    include_package_data=True,
    zip_safe=False,
)
