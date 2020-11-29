.. _install:

************
Installation
************

This chapter reviews a couple of different ways of installing PyGYRE
on your system.

Pre-requisites
==============

PyGYRE requires that the following Python packages are already present:

* `numpy <https://numpy.org/>`__
* `H5py <https://www.h5py.org/>`__
* `Astropy <https://www.astropy.org/>`__

If you opt to install from PyPI (below), then these pre-requisites
should be taken care of automatically.


Installing from PyPI
====================

To install PyGYRE from the `Python Package Index (PyPI)
<https://pypi.org/>`__, use the :command:`pip` command:

.. prompt:: bash

   pip install pygyre

If PyGYRE is already installed, you can upgrade to a more-recent
version via

.. prompt:: bash

   pip install --upgrade pygyre

Installing from Source
======================

To install PyGYRE from source, download the `source code
<github-tarball_>`__ and unpack it from the command line using the
:command:`tar` utility:

.. prompt:: bash
   :substitutions:

   tar xf pygyre-|release|.tar.gz

Then, change into the :file:`source` subdirectory of the newly created
directory and run the :command:`setup.py:` script:

.. prompt:: bash
   :substitutions:

   cd pygyre-|release|/source
   python setup.py install
   


	 



