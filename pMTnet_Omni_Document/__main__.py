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

from pMTnet_Omni_Document import __version__, __author__
from pMTnet_Omni_Document.data_curation import read_file


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
parser.add_argument("--output_folder_path",
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

    _, _ = read_file(file_path=file_path, 
                    background_tcrs_dir=background_tcrs_dir,
                    mhc_path=mhc_path,
                    save_results=True,
                    output_folder_path=cmdargs.output_folder_path,
                    sep=sep)
    
    sys.exit(0)


if __name__ == "__main__":
    cmdargs = parser.parse_args()
    main(cmdargs=cmdargs)
