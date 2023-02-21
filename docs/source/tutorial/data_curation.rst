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
        
        from pMTnet_Omni_Document.data_curation import read_file, encode_mhc_seq

Read the File 
--------------------
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

To perform data curation, simply call 

.. code-block:: python 

    df = read_file(file_path='file/path/here.csv',
                  background_tcrs_dir="./validation_data/",
                  mhc_path="./validation_data/valid_mhc.txt",
                  save_df=True,
                  output_file_path='output/file/path/here.csv',
                  sep=',')

You can inspect the values in ``df`` as well as all the 
``.csv`` files save to the path you specified with a few
modifications on the file names. 

Encode MHC sequences 
--------------------------
When you inspect the ``df``, if you see that some values 
in ``mhca_use_seq`` and/or ``mhcb_use_seq`` are ``True``, 
this means that for those records, we will use the sequences 
as either the MHCs are missing or we can not find their records
in our reference data. 

In this case, we will invoke the ESM2 algorithm to encode these 
sequences, produce a dictionary whose keys are the sequences 
and values are the embeddings, and save the dictionary 
as a ``.json`` file. 

To encode the MHC sequences, simply call 

.. code-block:: python 

    encode_mhc_seq(df=df, 
                   output_path='output/path.json')

.. warning:: 
    When uploading your dataset or the curated version produced 
    by ``read_file`` function (recommended)
    to `DBAI <http://lce-test.biohpc.swmed.edu/pmtnet>`_,
    make sure to upload the ``.json`` file as well. Otherwise, 
    these records will be dropped. 
