import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PALS-pathway",
    version="1.0.2",
    author="Joe Wandy",
    author_email="joe.wandy@glasgow.ac.uk",
    description="Pathway-level Analysis of Metabolites Expression Data",
    url="https://github.com/glasgowcompbio/PALS",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)