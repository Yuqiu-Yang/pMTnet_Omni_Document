Data Curation
===================
So I see you want to get into the nitty gritty details of our package. 
While our :doc:`../api_reference/index`
provides all the details of the functions defined in the package, 
it could be a bit overwhelming. In this relatively detailed tutorial, we will 
provide an in-depth exploration of how the package works. 

.. note:: 
    We will assume that the user has already imported 
    necessary functions using the command:

    .. code-block:: python
        
        from pMTnet_Omni_Document.data_curation import read_file


To perform data curation, the ``read_file`` function carries out 
the following steps sequentially: 

Check column names 
    Make sure all necessary columns are present 

Check species 
    Make sure values in ``tcr_species`` and ``pmhc_species``
    are ``human`` or ``mouse``

Check VA and Vb 
    If ``va`` and ``vb`` are not missing, we will look 
    up their corresponding sequences using our 
    reference data at ``./validation_data`` even if 
    users supplied the sequences. If we can not find 
    the gene/allele names, the records will be dropped.
    
    A separate dataframe ``./filename_curated_invalid_v.csv``
    will contain dropped records. 

    Otherwise, we will use the sequence information 
    with minimal curation 

Infer MHC information 
    Based on information provided in ``pmhc_species``, we 
    will infer the information on alpha and beta chains as 
    well as MHC classes 

Check MHC 
    If ``mhc`` is not missing, we will first look up its
    information in our reference data. If we can not 
    found the MHC or the MHC is missing,
    we will use the sequence information.

    Two separate dataframes 
    ``./filename_curated_mhc_alpha_dropped.csv`` and 
    ``./filename_curated_mhc_beta_dropped.csv``
    will contain dropped records. 
    
Check peptide 
    Sequences containing more than 30 AAs will be dropped. 

    A separate dataframe ``./filename_curated_antigen_dropped.csv``
    will contain dropped records. 

Check columns with amino acids
    Unknown AAs will be replaced with ``_``

Encode MHC sequences 
    When the ``df`` contains some MHCs that are NOT in our reference data or 
    are missing, we will inzoke the ESM2 algorithm to encode these 
    sequences, produce a dictionary whose keys are the sequences 
    and values are the embeddings, and save the dictionary 
    as a ``.json`` file. 

To perform data curation, simply call 

.. code-block:: python 

    df = read_file(file_path='file/path/here.csv',
                  background_tcrs_dir="./validation_data/",
                  mhc_path="./validation_data/valid_mhc.txt",
                  save_df=True,
                  output_folder_path='output/folder/path/',
                  sep=',')

You can inspect the values in ``df`` as well as all the 
``.csv`` files save to the path you specified with a few
modifications on the file names. 


.. warning:: 
    When uploading your dataset or the curated version produced 
    by ``read_file`` function (recommended)
    to `DBAI <http://lce-test.biohpc.swmed.edu/pmtnet>`_,
    make sure to upload the ``.json`` file as well. Otherwise, 
    these records will be dropped. 

.. warning:: 
    The current version of **pMTnet_Omni** will **NOT** perform 
    data curation. It will read in the data as is. If the data 
    do not conform with its required format, the program WILL halt. 