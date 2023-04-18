Quick Start Guide
==================
So I hear your week has been hectic. 
No worries. In this tutorial, we will walk through the basic usage of 
**pMTnet Omni Document** with minimum configuration. 
If you are truly swamped, we recommend `our online tool <http://lce-test.biohpc.swmed.edu/pmtnet>`_.


.. note::
    We will use the validation data we published along with the package, whose location we 
    assume is at ``./validation_data``. 
    
    Make sure your data file which we assume is located at 
    ``./df.csv`` is structured somewhat like the following:

    .. image:: ./images/sample_df.png
        :width: 600
        :align: center

    For a more detailed instruction on the data format, please check out :doc:`input_format/index`. 


CLI (Command Line Interface)
--------------------------------
By using CLI, you only need one line of code. 

.. code:: bash 

    python -m pMTnet_Omni_Document --file_path ./df.csv --validation_data_path ./validation_data --output_folder_path ./


Interactive Python 
-------------------
Read the file 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python 

    # Import necessary functions
    from pMTnet_Omni_Document.data_curation import read_file 

    # Read the file 
    df, mhc_seq_dict = read_file(file_path="./df.csv",
                                background_tcrs_dir='./validation_data/',
                                mhc_path="./validation_data/valid_mhc.txt",
                                save_results=True,
                                sep=",")

``./df_curated.csv`` will contain all the curated data. You will 
also see some extra columns in this file. 

``./mhc_seq_dict.json`` is a *json* file of a dictionary.
The keys are various MHC sequences and the values are their corresponding 
ESM embeddings.

.. note:: 
    If you see that ``mhca_use_seq`` and/or ``mhcb_user_seq`` columns 
    all have ``False``, then the json file will simply contain an empty 
    dictionary. 











