"""
Setup configuration for North Atlantic SST Pattern Analysis
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="north_atlantic_sst",
    version="0.1.0",
    author="Quan Liu",
    author_email="",
    description="Tools for analyzing SST patterns in the North Atlantic",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/liuquan18/North_Atlantic_SST_pattern",
    packages=find_packages(where="."),
    package_dir={"": "."},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.24.0",
        "xarray>=2023.1.0",
        "scipy>=1.10.0",
        "pandas>=2.0.0",
        "netcdf4>=1.6.0",
        "matplotlib>=3.7.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
        ],
        "viz": [
            "cartopy>=0.21.0",
            "seaborn>=0.12.0",
        ],
    },
)
