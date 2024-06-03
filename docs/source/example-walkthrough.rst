.. _walkthrough:

*******************
Example Walkthrough
*******************

This chapter walks through a number of typical PyGYEE use-cases.

Reading GYRE Output
===================

PyGYRE's most important capability is reading the :ref:`summary
<gyre:summary-files>` and :ref:`detail <gyre:detail-files>` files
created by GYRE. Assuming you're in the working directory you created
for the GYRE :ref:`example walkthrough <gyre:walkthrough>`, this
illustrates how you can read and inspect the summary file using the
:func:`pygyre.read_output` function:

.. code::

   import pygyre
   import numpy

   # Read data from the summary file

   s = pygyre.read_output('summary.h5')

   # Print the table metadata and contents

   print(s.meta)
   print(s)

In this example, `s` is a :py:class:`astropy.table.Table` instance, and
therefore supports a rich set of data analysis and processing operations. For instance,

.. code::

   # Print only the rows for modes with harmonic degree l = 2

   print(s[s['l'] == 1])

   # Print the indices of the mode with the largest inertia

   i = numpy.argmax(s['E_norm'])

   print(s[('l','n_pg')][i])

The same :func:`pygyre.read_output` function can also be used to read
detail files:

.. code::

   # Read data from one of the detail files

   d = pygyre.read_output('detail.l1.n-5.h5')

   # Print the table metadata and contents

   print(d.meta)
   print(d)

Reading Stellar Models
======================

PyGYRE can read stellar models the :ref:`MESA
<gyre:mesa-file-format>`, :ref:`GSM <gyre:gsm-file-format>` and
:ref:`POLY <gyre:poly-file-format>` formats. Again assuming that
you're in the working directory for the GYRE :ref:`example walkthrough
<gyre:walkthrough>`, here's how you can read the model file for the
:math:`5\,{\rm M_{\odot}}` slowly pulsating B star:

.. code::

   import pygyre

   # Read data from the model file

   m = pygyre.read_model('spb.mesa')

   # Print the table metadata and contents

   print(m.meta)
   print(m)
