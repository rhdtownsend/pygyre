import setuptools
import subprocess
import re

with open("README.md", "r") as fh:
    long_description = fh.read()

version_str = subprocess.check_output(["git", "describe", "--tags", "--abbrev=0"], text=True)
version = re.sub(r'^v', r'', version_str)

setuptools.setup(
    name="pygyre",
    version=version,
    author="Rich Townsend",
    author_email="townsend@astro.wisc.edu",
    description="Python support for the GYRE stellar oscillation code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rhdtownsend/pygyre",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    install_requires=[
        "h5py",
        "astropy",
        "numpy"
    ],
    python_requires='>=3.6'
)
