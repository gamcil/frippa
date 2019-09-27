from setuptools import setup, find_packages

with open("README.md") as readme:
    long_description = readme.read()

setup(
    name="fRiPPa",
    author="Cameron Gilchrist",
    version="0.0.1",
    description="fungal Ribosomally synthesized and Post-translationally"
    " modified Peptide (RiPP) assayer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gamcil/frippa",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["frippa=frippa.main:main"]},
)
