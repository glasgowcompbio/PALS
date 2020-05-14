from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="PALS-pathway",
    version="1.3.0",
    author="Joe Wandy",
    author_email="joe.wandy@glasgow.ac.uk",
    description="Pathway-level Analysis of Metabolites Expression Data",
    url="https://github.com/glasgowcompbio/PALS",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=['pals/run.py'],
    python_requires='>=3.6',
    packages=find_packages(),
    package_data={
        'pals': [
            'data/*.json.zip',
            'data/reactome/*/*/*.json.zip'
        ],
    }
)