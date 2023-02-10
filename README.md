# pMTnet Omni: your one-stop TCR-pMHC affinity prediction algorithm :microscope:

![Logo](/assets/logo.png)

![forthebadge](/assets/tcr-pmhc.svg)
![forthebadge](/assets/deep-learning.svg)

[![Documentation Status](https://readthedocs.org/projects/pmtnet-omni-document/badge/?version=latest)](https://pmtnet-omni-document.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni/branch/main/graph/badge.svg?token=L59TPMM3VN)](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni)

<b>pMTnet Omni</b> is a deep learning algorithm for affinity prediction based on TCR Va, Vb, CDR3a, CDR3b sequences, peptide sequence, and MHC allele types. The predictions can be made for human and mouse alleles, and for both CD8 T cells/MHC class I and CD4 T cells/MHC class II.

Please refer to our paper for more details: [pMTnet Omni paper link here](www.google.com)

We host the online tool on [DBAI](http://lce-test.biohpc.swmed.edu/pmtnet), where you can find all the members of the pMTnet 
family, including [pMTnet V1](https://github.com/tianshilu/pMTnet). 

We have also built a detailed [online documentation](https://pmtnet-omni-document.readthedocs.io/en/latest/) where we guide you step-by-step on how to format your data so it can be accpted by our algorithm.

## Model Overview 
![Model Overview](/docs/source/images/overview.png)


## Quick Start Guide 
1. Prepare your dataset so that it looks somewhat like the following:
![Sample df](/docs/source/images/sample_df.png)
Along with the main program, we also published 5 datasets under the `./validation_data` folder. Feel free 
to use those datasets to check if you TCR namings, Amino Acid sequences, and MHC namings conform with our 
standard.

**_NOTE:_** If the corresponding AA sequences are NOT provided AND the TCR names (same for MHC) can NOT be found in these datasets, the record WILL be dropped. 

2. Go to [our website](http://lce-test.biohpc.swmed.edu/pmtnet) and upload your data

3. An example output would look like this:
![Sample output](/docs/source/images/sample_output.png)

For a more in-depth explanation on input format, check out our [online documentation](https://pmtnet-omni-document.readthedocs.io/en/latest/). 

## CITATION HERE 

