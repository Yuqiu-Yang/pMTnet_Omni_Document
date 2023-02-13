"""The Command-Line Interface (CLI) of pMTnet Omni Document

The CLI of pMTnet Omni can be accessed via ``python -m pMTnet_Omni_Document``.

:Example:

    Get help:
    
    .. code-block:: bash

        python -m pMTnet_Omni -h
    
    Check version and authors:
    
    .. code-block:: bash
    
        python -m pMTnet_Omni_Document --version 
        python -m pMTnet_Omni_Document --author

**Pros**
    * Shallow learning curve
    
**Cons**
    * Fewer configurations 
    * No inspections of intermediate results
"""


import os
import sys
import argparse
import pandas as pd

from pMTnet_Omni_Document import __version__, __author__
from pMTnet_Omni_Document.data_curation import check_column_names,\
    check_species,\
    check_va_vb,\
    infer_mhc_info,\
    check_mhc,\
    encode_mhc_seq,\
    check_peptide,\
    check_amino_acids_columns


parser = argparse.ArgumentParser(description="pMTnet Omni")

parser.add_argument("--version", action="version",
                    version=__version__, help="Display the version of the software")
parser.add_argument("--author", action="version", version=__author__,
                    help="Check the author list of the algorithm")

# User input
parser.add_argument("--file_path", help="The path to the user data file")
parser.add_argument("--validation_data_path",
                    help="The path to the validation data files")

# User output
parser.add_argument("--output_file_path",
                    help="The file path to the output file")

# Maybe consider
# config.yaml


def main(cmdargs: argparse.Namespace):
    """The main method for pMTnet Omni

    Parameters:
    ----------
    cmdargs: argparse.Namespace
        The command line argments and flags 
    """

    file_path = cmdargs.file_path

    file_extention = os.path.splitext(file_path)[1]
    if file_extention in [".txt", ".tsv"]:
        sep = "\t"
    elif file_extention == ".csv":
        sep = ","
    else:
        raise Exception("Only .txt, .tsv, and .csv are accepted.")

    background_tcrs_dir = cmdargs.validation_data_path
    mhc_path = os.path.join(cmdargs.validation_data_path, "valid_mhc.txt")

    output_file_name = os.path.splitext(cmdargs.output_file_path)[0]
    if os.path.splitext(cmdargs.output_file_path)[1] != ".csv":
        raise Exception(
            "Currently, only .csv files can be used as output files")

    # Read user data file
    df = pd.read_csv(file_path, sep=sep).fillna('')
    print("Number of rows in raw dataset: " + str(df.shape[0]))
    # Check column names
    df = check_column_names(df=df)
    # Check species
    df = check_species(df=df)
    # Check VA VB
    df, invalid_v_df = check_va_vb(df=df, background_tcrs_dir=background_tcrs_dir)
    # Check MHC
    df = infer_mhc_info(df=df)
    df, df_mhc_alpha_dropped, df_mhc_beta_dropped = check_mhc(df=df, mhc_path=mhc_path)
    # Check peptide
    df, df_antigen_dropped = check_peptide(df=df)
    # Check aa sequences
    df = check_amino_acids_columns(df=df)
    
    output_path = output_file_name+"_mhc_seq_dict.pickle"
    encode_mhc_seq(df=df, output_path=output_path)
    df.to_csv(cmdargs.output_file_path, sep=',', index=False)
    invalid_v_df.to_csv(output_file_name+"_invalid_v.csv", sep=',', index=False)
    df_mhc_alpha_dropped.to_csv(output_file_name+"_mhc_alpha_dropped.csv", sep=',', index=False)
    df_mhc_beta_dropped.to_csv(output_file_name+"_mhc_beta_dropped.csv", sep=',', index=False)
    df_antigen_dropped.to_csv(output_file_name+"_antigen_dropped.csv", sep=',', index=False)
    
    sys.exit(0)


if __name__ == "__main__":
    cmdargs = parser.parse_args()
    main(cmdargs=cmdargs)
