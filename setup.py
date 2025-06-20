#!/usr/bin/env python3
"""
Setup script for ViRAnPy - Viral Metagenomic Analysis Pipeline
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "ViRAnPy - Viral Metagenomic Analysis Pipeline"

# Read requirements
def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(requirements_path):
        with open(requirements_path, 'r') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return []

setup(
    name="viranpy",
    version="0.1.0",
    author="Naveen Duhan",
    author_email="naveen.duhan@outlook.com",
    description="Viral Metagenomic Analysis Pipeline - Comprehensive viral genome assembly, annotation, and analysis",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/naveenduhan/viranpy",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
    },
    entry_points={
        "console_scripts": [
            "viranpy=viranpy.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "viranpy": ["databases/*", "databases/*/*"],
    },
    zip_safe=False,
) 