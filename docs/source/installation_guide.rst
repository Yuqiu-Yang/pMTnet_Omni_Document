Installation Guide 
===================

Dependencies
-------------
The dependencies of pMTnet Omni Document is fairly standard for a deep learning-based application

.. list-table:: Dependencies
  :widths: 50 50
  :align: center
  :header-rows: 1

  * - Package 
    - Version 
  * - python
    - ``>=3.9``  
  * - numpy
    - ``==1.22.4``
  * - pandas
    - ``==1.5.2``
  * - tqdm 
    - ``==4.64.1``
  * - torch
    - ``==1.13.1``
  * - fair-esm
    - ``==2.0.0``

Environment Setup 
------------------
.. _installation guide:

Of course, you can create your own environment and install those dependencies manually. 
However, we have made this simple for you.

* Go to `our website <http://lce-test.biohpc.swmed.edu/pmtnet>`_ and download the environment files. Or, run 
  
  .. code:: bash 
    
    git clone https://github.com/Yuqiu-Yang/pMTnet_Omni_Document.git
* CD to the directory that contains both `pMTnet_Omni_Document_env.yml` and `requirements.txt` files
* Run the followig  
  
  .. code:: bash 
    
    conda env create -f pMTnet_Omni_Document_env.yml

This will create a conda environment *pMTnet_Omni_Document*

Package Installation
------------------------

.. code:: bash

   conda activate pMTnet_Omni_Document
   pip install pMTnet_Omni_Document
  
To quickly test if it has been installed:

.. code:: bash

  python -m pMTnet_Omni_Document --version 


