import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='nnk',
    version='0.0.1dev',
    #packages=['nnk',],
    description='Tools for counting "NNK" variants from a yeast display library and downstream analysis',
    author="CollinJ0",
    url= 'https://www.github.com/CollinJ0/nnk',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    scripts=['bin/nnk'],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
