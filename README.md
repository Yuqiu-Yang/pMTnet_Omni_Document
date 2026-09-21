# pMTnet Omni: pan-MHC and cross-Species Prediction of T Cell Receptor-Antigen Binding :microscope:

![Logo](/assets/pmtnet_logo.png)

![forthebadge](/assets/tcr-pmhc.svg)
![forthebadge](/assets/deep-learning.svg)

| Package | Documentation | Code Coverage |
| --- | --- | --- |
| pMTnet Omni | | [![codecov](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni/branch/main/graph/badge.svg?token=L59TPMM3VN)](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni) |
| pMTnet Omni Document | [![Documentation Status](https://readthedocs.org/projects/pmtnet-omni-document/badge/?version=latest)](https://pmtnet-omni-document.readthedocs.io/en/latest/?badge=latest) | [![codecov](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni_Document/branch/main/graph/badge.svg?token=BR1vyICN3q)](https://codecov.io/gh/Yuqiu-Yang/pMTnet_Omni_Document) |

<b>pMTnet Omni</b> is a deep learning algorithm for affinity prediction based on TCR Va, Vb, CDR3a, CDR3b sequences, peptide sequence, and MHC allele types. The predictions can be made for human and mouse alleles, and for both CD8 T cells/MHC class I and CD4 T cells/MHC class II.

Please refer to our paper for more details: [pMTnet Omni paper link here](https://doi.org/10.1038/s41467-026-73396-3)

We have also built a detailed [online documentation](https://pmtnet-omni-document.readthedocs.io/en/latest/) where we guide you step-by-step on how to format your data so it can be accpted by our algorithm.

**_NOTE:_** This is the documentation for the data curation supporting tool for <b>pMTnet Omni</b>. Use this BEFORE you upload your dataset to [DBAI](https://dbai.biohpc.swmed.edu/pmtnet/).

## Model Overview 
![Model Overview](/docs/source/images/overview.png)

## DBAI :computer:
To try out pMTnet Omni, we recommend our online tool hosted on [DBAI](https://dbai.biohpc.swmed.edu/pmtnet/), where you can find all the members of the pMTnet 
family, including [pMTnet V1](https://github.com/tianshilu/pMTnet). 

**_NOTE:_** Just upload the data that conforms with our input requirements to our server and it will curate the data and crunch the numbers for you. No need to use pMTnet_Omni_Document for curation. 

<details>
<summary>DIY </summary>
## DIY :muscle: 
### Dependencies 
- numpy==1.22.4
- pandas==1.5.2
- tqdm==4.64.1
- torch==1.13.1
- fair-esm==2.0.0

### Enviroment Setup
```shell
conda env create -f pMTnet_Omni_Document_env.yml
```

### Installation 
```shell
conda activate pMTnet_Omni_Document
pip install pMTnet_Omni_Document
```

</details>

### Quick Start Guide 
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

<details>
<summary> DIY </summary>

2. Say your dataset is under <i>./df.csv</i>. In your terminal, run 
```shell
conda activate pMTnet_Omni_Document

python -m pMTnet_Omni_Document --file_path ./df.csv --output_folder_path ./
```

3. Go to [our website](https://dbai.biohpc.swmed.edu/pmtnet/) and upload your data including the `.json` file. 

4. An example output would look like this:
![Sample output](/docs/source/images/sample_output.png)

</details>

For a more in-depth explanation on input format, check out our [online documentation](https://pmtnet-omni-document.readthedocs.io/en/latest/). 


### CITATION
We have uploaded our article to bioRxiv. To cite
```
@ARTICLE{Han2026-wo,
  title     = "Deciphering small sequence differences in {T} cell
               receptor-antigen pairing",
  author    = "Han, Yi and Yang, Yuqiu and Zhu, James and Fattah, Farjana J and
               von Itzstein, Mitchell S and Zhang, Minying and Bermack, Casey
               and Jiang, Peixin and Singh, Shailbala and Tian, Yanhua and Hu,
               Yifei and Deng, Yafang and Kang, Xiongbin and Yang, Donghan M
               and Liu, Jialiang and Xue, Yaming and Liang, Chaoying and Raman,
               Indu and Zhu, Chengsong and Xiao, Olivia and Dowell, Jonathan E
               and Homsi, Jade and Rashdan, Sawsan and Pan, Ke and Yang,
               Shengjie and Gwin, Mary E and Hsiehchen, David and
               Gloria-McCutchen, Yvonne and Wu, Fangjiang and Heymach, John V
               and Gibbons, Don and Huang, Junzhou and Cheng, Chao and Zhang,
               Jianjun and Yee, Cassian and Reuben, Alexandre and Gerber, David
               E and Wang, Tao",
  abstract  = "T cells have important functions in development and disease
               processes through T cell receptor (TCR)-dependent activities.
               Many tools were developed to predict the binding between TCRs
               and antigens. However, one of the uncertainties is whether such
               tools can decipher how small changes in the TCRs or antigenic
               peptides contribute to binding. We develop a deep learning
               model, pMTnet-omni, which not only predicts the binding vs.
               non-binding of TCRs towards pMHCs, but also distinguishes the
               stronger vs. weaker binding of TCRs similar in sequence. We
               leverage this capability to interpret the biological rules that
               govern TCR-antigen pairing. This also enables pMTnet-omni to
               accurately predict variant TCRs with desired stronger or weaker
               binding to the antigen, in conjunction with a Lab-in-the-Loop
               (LiL) mechanism. We show that pMTnet-omni can also predict
               binding of TCRs towards similar pMHCs. Overall, we provide a
               flexible toolkit for research and translational applications
               involving antigens and TCRs.",
  journal   = "Nat. Commun.",
  publisher = "Springer Science and Business Media LLC",
  volume    =  17,
  number    =  1,
  month     =  jun,
  year      =  2026,
  copyright = "https://creativecommons.org/licenses/by-nc-nd/4.0",
  language  = "en"
}

```



