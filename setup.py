from setuptools import setup, find_packages

setup(
    name="dna_dev_kit",
    version="0.2.0",
    description="DNA sequence analysis, restriction digest, primer design, codon optimisation, and Golden Gate assembly.",
    author="dna-dev-kit contributors",
    packages=find_packages(exclude=["tests", "tests.*"]),
    python_requires=">=3.9",
    install_requires=["numpy>=1.24"],
    extras_require={"dev": ["pytest>=7.0"]},
)
