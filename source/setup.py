import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pygyre",
    version="1.1",
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
