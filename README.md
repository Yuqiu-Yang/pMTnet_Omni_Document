# pMTnet Omni: your one-stop TCR-pMHC affinity prediction algorithm :microscope:

![Logo](/assets/logo.png)

![forthebadge](/assets/tcr-pmhc.svg)
![forthebadge](/assets/deep-learning.svg)

| Package | Documentation | Code Coverage |
| --- | --- | --- |
| pMTnet Omni | | [![codecov](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni/branch/main/graph/badge.svg?token=L59TPMM3VN)](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni) |
| pMTnet Omni Document | [![Documentation Status](https://readthedocs.org/projects/pmtnet-omni-document/badge/?version=latest)](https://pmtnet-omni-document.readthedocs.io/en/latest/?badge=latest) | [![codecov](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni_Document/branch/main/graph/badge.svg?token=BR1vyICN3q)](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni_Document) |

<b>pMTnet Omni</b> is a deep learning algorithm for affinity prediction based on TCR Va, Vb, CDR3a, CDR3b sequences, peptide sequence, and MHC allele types. The predictions can be made for human and mouse alleles, and for both CD8 T cells/MHC class I and CD4 T cells/MHC class II.

Please refer to our paper for more details: [pMTnet Omni paper link here](www.google.com)

We host the online tool on [DBAI](http://lce-test.biohpc.swmed.edu/pmtnet), where you can find all the members of the pMTnet 
family, including [pMTnet V1](https://github.com/tianshilu/pMTnet). 

We have also built a detailed [online documentation](https://pmtnet-omni-document.readthedocs.io/en/latest/) where we guide you step-by-step on how to format your data so it can be accpted by our algorithm.

**_NOTE:_** This is the documentation for the data curation supporting tool for <b>pMTnet Omni</b>. Use this BEFORE you upload your dataset to [DBAI](http://lce-test.biohpc.swmed.edu/pmtnet).

## Model Overview 
![Model Overview](/docs/source/images/overview.png)

## Dependencies 
- numpy==1.22.4
- pandas==1.5.2
- tqdm==4.64.1
- torch==1.13.1
- fair-esm==2.0.0

## Enviroment Setup
```shell
conda env create -f pMTnet_Omni_Document_env.yml
```

## Installation 
```shell
conda activate pMTnet_Omni_Document
pip install pMTnet_Omni_Document
```

## Quick Start Guide 
1. Prepare your dataset so that it looks somewhat like the following:
![Sample df](/docs/source/images/sample_df.png)
Along with the main program, we also published 5 datasets under the `./validation_data` folder. Feel free 
to use those datasets to check if you TCR namings, Amino Acid sequences, and MHC namings conform with our 
standard.

**_NOTE:_** When both TCR names (resp. MHC) and the 
TCR sequences (resp. MHC sequences) are provided, we 
will *disregard the sequences*. If the names can NOT be 
found in our reference database, the record WILL be 
dropped.

**_NOTE:_** On the other hand, if the names are NOT provided, we will use the sequences with minimal curation. 

2. Say your dataset is under <i>./df.csv</i>. In your terminal, run 
```shell
conda activate pMTnet_Omni_Document

python -m pMTnet_Omni_Document --file_path ./df.csv --validation_data_path ./validation_data --output_file_path ./df_result.csv
```

3. Go to [our website](http://lce-test.biohpc.swmed.edu/pmtnet) and upload your data including the `.pickle` file. 

4. An example output would look like this:
![Sample output](/docs/source/images/sample_output.png)

For a more in-depth explanation on input format, check out our [online documentation](https://pmtnet-omni-document.readthedocs.io/en/latest/). 

## CITATION HERE 

