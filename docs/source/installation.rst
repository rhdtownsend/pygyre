.. _install:

************
Installation
************

This chapter reviews a couple of different ways of installing PyGYRE
on your system.

Pre-requisites
==============

PyGYRE requires Python 3.9 or later, together with the following
Python packages:

* `numpy <https://numpy.org/>`__ (version 1.21 or later; NumPy 2 is supported)
* `H5py <https://www.h5py.org/>`__ (version 3.9 or later)
* `Astropy <https://www.astropy.org/>`__ (version 5.3.2 or later)

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

.. code-block:: console
   :substitutions:

   $ tar xf pygyre-|release|.tar.gz

Then, change into the newly created directory and install the package
with :command:`pip`:

.. code-block:: console
   :substitutions:

   $ cd pygyre-|release|
   $ pip install .
