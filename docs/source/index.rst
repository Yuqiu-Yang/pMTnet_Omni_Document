.. pMTnet_Omni documentation master file, created by
   sphinx-quickstart on Wed Dec 14 11:36:34 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ../../assets/logo.png
   :width: 600

**pMTnet Omni**: your one-stop TCR-pMHC affinity prediction algorithm
===========================================================================

**pMTnet Omni** is a deep learning algorithm for affinity prediction based on 
TCR Va, Vb, CDR3a, CDR3b sequences, peptide sequence, and MHC allele types.
The predictions can be made for human and mouse alleles, 
and for both CD8 T cells/MHC class I and CD4 T cells/MHC class II.

Here is a quick overview of the structure of the model:

.. image:: ./images/overview.png
   :width: 600
   :align: center

For a more detailed exploration of our model, please refer to `our paper <https://www.google.com>`_:

*pMTnet Omni Paper citation here*

The online tool are hosted on `DBAI <http://lce-test.biohpc.swmed.edu/pmtnet>`_ where you can upload
your own dataset and we will crunch the numbers for you. 

In this document, we will mainly focus on two aspects:

Input Format
   A series of detailed explanations as well as functions will be provided to help you organize your dataset so that 
   the input can be correctly recognized by **pMTnet Omni**.

Input Parameters 
   **pMTnet Omni** can be configured in a various ways to suit your own need. We will walk through 
   the parameters you can use to alter the behaviors of our algorithm.


User Guide / Tutorial
=====================
.. toctree::
   :maxdepth: 1
   :caption: Get Started

   installation_guide
   quick_start
   
.. toctree::
   :maxdepth: 3
   :caption: Input Format

   input_format/index

.. toctree::
   :maxdepth: 1
   :caption: Detailed Tutorials

   tutorial/data_curation
   tutorial/details

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api_reference/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
