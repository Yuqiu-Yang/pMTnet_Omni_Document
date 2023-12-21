Quick Start Guide
==================
So I hear your week has been hectic. 
No worries. In this tutorial, we will walk through the basic usage of 
**pMTnet Omni Document** with minimum configuration. 
If you are truly swamped, we recommend `our online tool <https://dbai.biohpc.swmed.edu/pmtnet/>`_.

.. note::
    
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

    python -m pMTnet_Omni_Document --file_path ./df.csv --output_folder_path ./


Interactive Python 
-------------------
Read the file 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python 

    # Import necessary functions
    from pMTnet_Omni_Document.data_curation import read_file 

    # Read the file 
    df, mhc_seq_dict = read_file(file_path="./df.csv",
                                save_results=True,
                                sep=",")

In the output, you will see two files:

``./df_curated.csv`` will contain all the curated data. You will 
also see some extra columns in this file. 

.. note:: 
    
    If you see that ``mhca_use_seq`` and/or ``mhcb_user_seq`` columns 
    all have ``False``, then the json file will simply contain an empty 
    dictionary. 

``./mhc_seq_dict.json`` is a *json* file of a dictionary.
The keys are various MHC sequences and the values are their corresponding 
ESM embeddings.

.. note:: 

    If you want to reproduce the training and validation results 
    in our paper, you can download 
    `our training data <https://365utsouthwestern-my.sharepoint.com/:x:/g/personal/yuqiu_yang_utsouthwestern_edu/EYBVZgSOuq9HpYNKrb3c5jIB_87GlX5prJ6hNQdaB77ltw?e=36ZmW0>`_ and 
    `our validation data <https://365utsouthwestern-my.sharepoint.com/:x:/g/personal/yuqiu_yang_utsouthwestern_edu/Eb4tEGiZS4tDtRUvlu_IoKYBf-zm_aojgo0tdVbTl2Au7Q?e=1ezESp>`_.
