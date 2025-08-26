"""Setup script for SeqMat"""
from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="seqmat",
    version="0.1.29",
    author="Nicolas Lynn Vila",
    author_email="nicolasalynn@gmail.com",
    description="Lightning-fast gene manipulation and analysis library.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nicolasalynn/seqmat",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.10",
    install_requires=[
        "numpy>=1.20.0,<2.0",
        "pandas>=1.3.0,<2.1.0",
        "pysam>=0.19.0",
        "requests>=2.26.0",
        "tqdm>=4.62.0",
        "biopython>=1.79",
        "gtfparse>=1.2.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.910",
        ],
    },
    entry_points={
        "console_scripts": [
            "seqmat=seqmat.cli:main",
            "seqmat-setup=seqmat.cli:main",  # Legacy alias
        ],
    },
    package_data={
        "seqmat": ["*.json"],
    },
    include_package_data=True,
)