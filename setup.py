import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='NNK',
    version='0.0.1dev',
    #packages=['NNK',],
    description='Tools for couting and analyzing NNK based deep mutational scanning libraries',
    author="CollinJ0",
    url= 'https://www.github.com/CollinJ0/NNK',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    scripts=['bin/NNK-counter', 'bin/NNK-analyzer'],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
