pMTnet Omni Tutorial
=================================
Although **pMTnet Omni** works well with its default setting,
we do provide the users with various parameters so that 
our tool can be configured to cater to your needs. 

In this section, we will provide a relatively detailed 
explanation on important parameters that you can tweak.

Step 1: Instantiation
-----------------------------------------------------
**pMTnet Omni** starts by instantiate a ``pMTnet_Omni_class`` object.
There are only two parameters needed to accomlish this. 

.. list-table:: Instantiation Parameters 
    :align: center
    :widths: 10 10 20 10 50 
    :header-rows: 1

    * - Parameter 
      - Input Format 
      - Acceptable Inputs 
      - Default
      - Note
    * - model_device 
      - ``str``
      - "cpu", "gpu", None
      - None
      - If ``None``, choose based on gpu availability
    * - seed 
      - ``int``
      - Any integer, None
      - None 
      - If ``None``, no seed is used

Step 2: Data Importing 
------------------------
During our internal testing, we found that while a mediocre GPU is capable 
of handling several to hundreds of pairs, for dataset containing thousands of 
TCR-pMHC pairs, we usually will run out of memory. Partitioning the original 
data frame is a simple yet effective remedy. Therefore, we 
provide a parameter ``partition_size`` for this exact purpose. 

.. image:: ../images/partitions.png
    :width: 600
    :align: center

Suppose the user's data contains 5000 TCR-pMHC pairs. 
By setting ``partition_size=2000`` we will partition the 
original data frame into three partitions with 
2,000, 2,000, and 1,000 pairs, respectively.

.. list-table:: Data Imporing Parameters 
    :align: center
    :widths: 10 10 20 10 50 
    :header-rows: 1

    * - Parameter 
      - Input Format 
      - Acceptable Inputs 
      - Default
      - Note
    * - partition_size 
      - ``int``
      - Any integer
      - 500
      - If ``model_device`` is "cpu", it's not really used

Step 3: TCR-pMHC Affinity Prediction
---------------------------------------
There are quite a lot you can configure at this stage. 

.. list-table:: Affinity Prediction Parameters 
    :align: center
    :widths: 10 10 20 10 50 
    :header-rows: 1

    * - Parameter 
      - Input Format 
      - Acceptable Inputs 
      - Default
      - Note
    * - compute_percentile_rank 
      - ``bool``
      - True, False 
      - False 
      - If False, only the raw affinity scores will be reported. 
    * - rank_threshold 
      - ``float``
      - Any number between 0 and 1
      - 0.03
      - The rank percentile threshold greater than which further verification is NOT conducted 
    * - B
      - ``int``
      - Any positive integer
      - 1
      - Number of trials 
    * - check_size 
      - ``list``
      - List of positive integers 
      - [1000, 10000, 100000]
      - Explanations provided below 
    * - load_size
      - ``int``
      - Any positive integer
      - 1,000,000
      - Explanations provided below 
    * - minibatch_size
      - ``int``
      - Any positive integer
      - 50,000
      - Explanations provided below 

In general, to get the final percentile ranks for each one of the 
TCR-pMHC pairs you provided in the dataframe, two steps are involved:
**getting the raw affinity scores** and 
(you guessed it) **computing the percentile ranks**. 

Raw Affinity Scores 
~~~~~~~~~~~~~~~~~~~~~~~
If the argument ``compute_percentile_rank`` is ``False``, the 
program will halt. And only the raw affinity scores will 
be reported. The rest of the parameters won't matter. 

On the other hand, if ``compute_percentile_rank=True``, 
the program will proceed to the next stage. 

Percentile Ranks 
~~~~~~~~~~~~~~~~~~~~~~~~
In this stage, we will compare each pair with the background TCRs. As 
we have millions of background TCRs, it would be time consuming to check 
each pair against the entire database. Hence, in our implementation, we borrowed 
an idea from the literature of clinical trials.

Each TCR-pMHC pair will undergo one or several "trials", each with a 
sequence of `checks`. The procedure is conceptually simple: 

Within a "trial", each TCR-pMHC pair will be first checked against a small 
subset of the background TCRs. If the predicted binding is strong enough, we sample 
a larger subset of the background TCRs and check the given pair against them. 
The process repeats until either the pair falls out of the top rank list or 
it has been validated against enough background TCRs, at which point, the 
algorithm reports the final rank. 

.. image:: ../images/prediction.png
    :width: 600
    :align: center

Two parameters ``load_size`` and ``minibatch_size`` could be somewhat confusing. But they 
are implemented to further speed up the prediction process. 

``load_size`` is implemented so that for each trial, only that many background TCRs will 
be potentially used. This is because the background TCRs datasets are relatively large, 
meaning that initializing the dataloaders will be time consuming. 

``minibatch_size`` is how many background TCRs the dataloader will sample within a check_size.
For example, is the current check_size is 2,000 and the minibatch_size is 1,000. Then 
the dataloader will first load 1,000 TCRs, compute the rank, load another 1,000 TCRs, and 
update the rank. This will speed up the process as directly load, say 1,000,000 TCRs will 
be slow. 

Sample Output
~~~~~~~~~~~~~~

.. image:: ../images/sample_output.png
    :width: 700
    :align: center

We also provide other supporting files for your to download. 

Parameter Summary
-----------------------

.. list-table:: Parameter Summary
    :align: center
    :widths: 10 10 20 10 50 
    :header-rows: 1

    * - Parameter 
      - Input Format 
      - Acceptable Inputs 
      - Default
      - Note
    * - model_device 
      - ``str``
      - "cpu", "gpu", None
      - None
      - If ``None``, choose based on gpu availability
    * - seed 
      - ``int``
      - Any integer, None
      - None 
      - If ``None``, no seed is used
    * - partition_size 
      - ``int``
      - Any integer
      - 500
      - If ``model_device`` is "cpu", it's not really used
    * - compute_percentile_rank 
      - ``bool``
      - True, False 
      - False 
      - If False, only the raw affinity scores will be reported. 
    * - rank_threshold 
      - ``float``
      - Any number between 0 and 1
      - 0.03
      - The rank percentile threshold greater than which further verification is NOT conducted 
    * - B
      - ``int``
      - Any positive integer
      - 1
      - Number of trials 
    * - check_size 
      - ``list``
      - List of positive integers 
      - [1000, 10000, 100000]
      - 
    * - load_size
      - ``int``
      - Any positive integer
      - 1,000,000
      - 
    * - minibatch_size
      - ``int``
      - Any positive integer
      - 50,000
      - 