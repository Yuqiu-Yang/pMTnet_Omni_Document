# Data IO
import os
import pandas as pd

# String operation
import re

# Progress bar
from tqdm import tqdm

# Typing
from typing import Tuple


def check_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Check if the column names are correct 

    The main purpose of this function is to make sure that the dataframe 
    provided by the users contains necessary columns so that 
    it can be used by subsequent functions 

    This function will NOT create mhca, mhcb, mhcaseq, mhcbseq.
    It will keep the original column mhc, mhcseq.

    Parameters
    ----------
    df: pd.DataFrame
        A pandas dataframe containing pairing data 

    Returns
    ---------
    pd.DataFrame:
        A pandas dataframe with corrected column names

    """
    # Define the column names to be used
    common_column_names = ["va", "cdr3a", "vaseq", "vb", "cdr3b", "vbseq",
                           "peptide", "mhc", "mhcseq",
                           "tcr_species", "pmhc_species"]
    # Get rid of spaces, convert names to lowercase, remove special characters
    print("Checking column names...\n")
    df.columns = df.columns.str.replace(' ', '').str.lower().str.replace(
        "_", "").str.replace(r'\W', '', regex=True)
    df_cols = df.columns.tolist()
    for i in range(len(df_cols)):
        # We get rid of white spaces from each column
        df[df_cols[i]] = df[df_cols[i]].str.replace(' ', '')
        df_cols[i] = re.sub(r'(?<=tcr).*(?=species)', "_", df_cols[i])
        df_cols[i] = re.sub(r'(?<=pmhc).*(?=species)', "_", df_cols[i])
    df.columns = df_cols
    ##################################################
    # We first make sure that all required columns are present
    # in the dataframe
    ##################################################
    # Get names not presented in the user dataframe
    user_df_set = set(df_cols)
    diff_set = set(common_column_names) - user_df_set

    error_messages = "".join(["Column " + name + " can not be found.\n"
                              for name in ["cdr3a", "cdr3b", "peptide",
                                           "tcr_species", "pmhc_species"] if name in diff_set])

    if error_messages != "":
        raise Exception(error_messages)
    ##################################################
    # mhc mhcseq
    ##################################################
    if ("mhc" in diff_set) and ("mhcseq" in diff_set):
        raise Exception("At least one of mhc and mhcseq need to exist")
    if ("mhc" in diff_set) and ("mhcseq" not in diff_set):
        # if column mhc is missing but we have mhcseq
        df['mhc'] = ''
    if ("mhc" not in diff_set) and ("mhcseq" in diff_set):
        # if column mhcseq is missing but we have mhc
        df['mhcseq'] = ''
    ##################################################
    # va vaseq
    ##################################################
    if ("va" in diff_set) and ("vaseq" in diff_set):
        raise Exception("At least one of va and vaseq need to exist")
    if ("va" in diff_set) and ("vaseq" not in diff_set):
        # if column va is missing but we have vaseq
        df['va'] = ''
    if ("va" not in diff_set) and ("vaseq" in diff_set):
        # if column vaseq is missing but we have va
        df['vaseq'] = ''
    ##################################################
    # vb vbseq
    ##################################################
    if ("vb" in diff_set) and ("vbseq" in diff_set):
        raise Exception("At least one of vb and vbseq need to exist")
    if ("vb" in diff_set) and ("vbseq" not in diff_set):
        # if column vb is missing but we have vbseq
        df['vb'] = ''
    if ("vb" not in diff_set) and ("vbseq" in diff_set):
        # if column vbseq is missing but we have vb
        df['vbseq'] = ''

    return df


def check_species(df: pd.DataFrame) -> pd.DataFrame:

    # We first perform some basic data curation to left some burden off the users
    df["tcr_species"] = df["tcr_species"].str.lower().str.replace(' ', '').str.replace(
        "_", "").str.replace(r'\W', '', regex=True)
    df["pmhc_species"] = df["pmhc_species"].str.lower().str.replace(' ', '').str.replace(
        "_", "").str.replace(r'\W', '', regex=True)
    # Then we make sure that the species are either human or mouse
    if not set(df['tcr_species'].values).issubset(set(['human', 'mouse'])):
        raise Exception("tcr_species have to be human or mouse")
    if not set(df['pmhc_species'].values).issubset(set(['human', 'mouse'])):
        raise Exception("pmhc_species have to be human or mouse")

    return df


def check_va_vb(df: pd.DataFrame,
                background_tcrs_dir: str = "./validation_data/") -> Tuple[pd.DataFrame, pd.DataFrame]:

    ##################################################
    # The basic logic is
    # If va vb are provided, then we will look up
    # the corresponding sequences in our reference
    # data REGARDLESS of the presence of the sequence
    # Otherwise, we simply use the sequence
    ##################################################

    human_alpha_tcrs = pd.read_csv(os.path.join(
        background_tcrs_dir, "human_alpha.txt"), sep="\t", header=0)[["va", "vaseq"]].drop_duplicates()
    human_beta_tcrs = pd.read_csv(os.path.join(
        background_tcrs_dir, "human_beta.txt"), sep="\t", header=0)[["vb", "vbseq"]].drop_duplicates()
    mouse_alpha_tcrs = pd.read_csv(os.path.join(
        background_tcrs_dir, "mouse_alpha.txt"), sep="\t", header=0)[["va", "vaseq"]].drop_duplicates()
    mouse_beta_tcrs = pd.read_csv(os.path.join(
        background_tcrs_dir, "mouse_beta.txt"), sep="\t", header=0)[["vb", "vbseq"]].drop_duplicates()
    # We rename the *seq* columns
    human_alpha_tcrs = human_alpha_tcrs.rename(columns={'vaseq': 'vaseq_ref'})
    human_beta_tcrs = human_beta_tcrs.rename(columns={'vbseq': 'vbseq_ref'})
    mouse_alpha_tcrs = mouse_alpha_tcrs.rename(columns={'vaseq': 'vaseq_ref'})
    mouse_beta_tcrs = mouse_beta_tcrs.rename(columns={'vbseq': 'vbseq_ref'})
    # We separate the dataframe based on species
    human_df = df[df['tcr_species'] == "human"].reset_index(drop=True)
    mouse_df = df[df['tcr_species'] == "mouse"].reset_index(drop=True)
    # Before we start, we curate the data a little
    human_df['va'] = human_df['va'].str.replace(
        ' ', '').str.replace('.', '-').str.replace(':', '*')
    mouse_df['va'] = mouse_df['va'].str.replace(
        ' ', '').str.replace('.', '-').str.replace(':', '*')
    human_df['vb'] = human_df['vb'].str.replace(
        ' ', '').str.replace('.', '-').str.replace(':', '*')
    mouse_df['vb'] = mouse_df['vb'].str.replace(
        ' ', '').str.replace('.', '-').str.replace(':', '*')

    # We will left join the two sets of dataframes on va/vb
    # Then depending on whether or not the seq is '' or not
    # we will choose if our reference should be used
    human_df = human_df.merge(human_alpha_tcrs, on="va", how='left')
    human_df = human_df.merge(human_beta_tcrs, on="vb", how='left')
    mouse_df = mouse_df.merge(mouse_alpha_tcrs, on="va", how='left')
    mouse_df = mouse_df.merge(mouse_beta_tcrs, on="vb", how='left')
    # .._df now should have columns va vb vaseq vbseq vaseq_ref vbseq_ref along with others
    for i in tqdm(range(human_df.shape[0])):
        va = human_df.at[i, 'va']
        vb = human_df.at[i, 'vb']
        # If va/vb can be found we use the reference
        # Otherwise, we use the user input
        if va != '':
            human_df.at[i, 'vaseq'] = human_df.at[i, 'vaseq_ref']
        if vb != '':
            human_df.at[i, 'vbseq'] = human_df.at[i, 'vbseq_ref']

    human_df = human_df.drop(columns=['vaseq_ref', 'vbseq_ref']).fillna('')
    invalid_human_df = human_df[(human_df['vaseq'] == '') |
                                (human_df['vbseq'] == '')].reset_index(drop=True)
    human_df = human_df[(human_df['vaseq'] != '') &
                        (human_df['vbseq'] != '')].reset_index(drop=True)

    for i in tqdm(range(mouse_df.shape[0])):
        va = mouse_df.at[i, 'va']
        vb = mouse_df.at[i, 'vb']
        # If va/vb can be found we use the reference
        # Otherwise, we use the user input
        if va != '':
            mouse_df.at[i, 'vaseq'] = mouse_df.at[i, 'vaseq_ref']
        if vb != '':
            mouse_df.at[i, 'vbseq'] = mouse_df.at[i, 'vbseq_ref']
    
    mouse_df = mouse_df.drop(columns=['vaseq_ref', 'vbseq_ref']).fillna('')
    invalid_mouse_df = mouse_df[(mouse_df['vaseq'] == '') |
                                (mouse_df['vbseq'] == '')].reset_index(drop=True)
    mouse_df = mouse_df[(mouse_df['vaseq'] != '') &
                        (mouse_df['vbseq'] != '')].reset_index(drop=True)

    df = pd.concat([human_df, mouse_df], axis=0,
                   ignore_index=True).reset_index(drop=True)
    invalid_df = pd.concat([invalid_human_df, invalid_mouse_df],
                           axis=0, ignore_index=True).reset_index(drop=True)
    return df, invalid_df


def infer_mhc_info(df: pd.DataFrame) -> pd.DataFrame:
    """Infer the MHC classes and create columns mhca and mhcb

    The input df should be the output of the check_column_names function.

    Parameters
    ---------
    df: pd.DataFrame
        A pandas dataframe containing mhc, mhcseq, and pmhc_species

    Returns
    ---------
    pd.DataFrame
        A column of a pandas dataframe with the inferred MHC classes, 
        MHCs on the alpha chain and the beta chain 
    """
    df = df.reset_index(drop=True)
    df['mhc_class'] = ''
    df['mhca'] = ''
    df['mhcb'] = ''
    df['mhcaseq'] = ''
    df['mhcbseq'] = ''
    # We clean up the input a little bit
    # We remove HLA and H2 which will be added backed if necessary
    # so as to make sure the format is coherent
    df['mhc'] = df['mhc'].str.replace(' ', '').str.replace(
        '-', '').str.replace('HLA', '').str.replace('H2', '')
    df['mhcseq'] = df['mhcseq'].str.replace(' ', '').str.upper()

    # Human class I starts with A, B, or C
    human_class_i = ("A", "B", "C")
    # Class II starts with DP, DQ, or DR
    human_class_ii = ("DP", "DQ", "DR")
    # Mouse class I starts with D, K, L, Q
    mouse_class_i = tuple(["H-2-" + i for i in ["D", "K", "L", "Q"]])
    # Class II starts with IA, IE
    mouse_class_ii = tuple(["H-2-" + i for i in ["IA", "IE"]])

    for i in range(df.shape[0]):
        mhc = df.at[i, "mhc"]
        mhcseq = df.at[i, "mhcseq"]
        pmhc_species = df.at[i, "pmhc_species"]
        if pmhc_species == "human":
            if mhc.startswith(human_class_i):
                # For human class I, mhca will be the
                # input mhc. mhcb will be human_microglobulin
                df.at[i, 'mhc_class'] = 'human class i'
                df.at[i, 'mhca'] = mhc
                df.at[i, 'mhcb'] = "human_microglobulin"
                df.at[i, 'mhcaseq'] = mhcseq
                df.at[i, 'mhcbseq'] = ''
            elif mhc.startswith(human_class_ii):
                # For human class II, usually alpha beta will both
                # be provided
                # However, for DR, its possible that only beta is provided
                df.at[i, 'mhc_class'] = 'human class ii'
                mhc_list = mhc.split('/')
                mhcseq_list = mhcseq.split('/')
                if len(mhcseq_list) == 1:
                    mhcseq_list.append('')
                if (len(mhc_list) == 1) and (mhc.startswith("DR")):
                    df.at[i, 'mhca'] = "DRA*01:01"
                    df.at[i, 'mhcb'] = mhc
                    df.at[i, 'mhcaseq'] = ''
                    df.at[i, 'mhcbseq'] = mhcseq
                else:
                    df.at[i, 'mhca'] = mhc_list[0]
                    df.at[i, 'mhcb'] = mhc_list[1]
                    df.at[i, 'mhcaseq'] = mhcseq_list[0]
                    df.at[i, 'mhcbseq'] = mhcseq_list[1]
            elif (mhc == '') and ('/' in mhcseq):
                # MHC is missing users have to provide
                # seq on both chains separated by /
                df.at[i, 'mhc_class'] = 'human'
                mhc_list = mhc.split('/')
                mhcseq_list = mhcseq.split('/')
                df.at[i, 'mhca'] = ''
                df.at[i, 'mhcb'] = ''
                df.at[i, 'mhcaseq'] = mhcseq_list[0]
                df.at[i, 'mhcbseq'] = mhcseq_list[1]
            else:
                raise ValueError(
                    "Missing sequence info or class of " + mhc + " can not be determined")
        elif pmhc_species == "mouse":
            if mhc != '':
                mhc = 'H-2-'+mhc
            if mhc.startswith(mouse_class_i):
                # For mouse Class I, it's on the alpha
                # beta will be mouse_microglobulin
                df.at[i, 'mhc_class'] = 'mouse class i'
                df.at[i, 'mhca'] = mhc
                df.at[i, 'mhcb'] = 'mouse_microglobulin'
                df.at[i, 'mhcaseq'] = mhcseq
                df.at[i, 'mhcbseq'] = ''
            elif mhc.startswith(mouse_class_ii):
                # For mouse Class II, alpha and beta
                # will share the same mhc
                mhcseq_list = mhcseq.split('/')
                if len(mhcseq_list) == 1:
                    mhcseq_list.append('')
                df.at[i, 'mhc_class'] = 'mouse class ii'
                df.at[i, 'mhca'] = mhc+"_alpha"
                df.at[i, 'mhcb'] = mhc+"_beta"
                df.at[i, 'mhcaseq'] = mhcseq_list[0]
                df.at[i, 'mhcbseq'] = mhcseq_list[1]
            elif (mhc == '') and ('/' in mhcseq):
                # MHC is missing users have to provide
                # seq on both chains separated by /
                df.at[i, 'mhc_class'] = 'mouse'
                mhc_list = mhc.split('/')
                mhcseq_list = mhcseq.split('/')
                df.at[i, 'mhca'] = ''
                df.at[i, 'mhcb'] = ''
                df.at[i, 'mhcaseq'] = mhcseq_list[0]
                df.at[i, 'mhcbseq'] = mhcseq_list[1]
            else:
                raise ValueError(
                    "Missing sequence info or class of " + mhc + " can not be determined")
        else:
            raise ValueError(pmhc_species + " is not human or mouse")

    return df


def check_mhc(df: pd.DataFrame,
              mhc_path: str = "./validation_data/valid_mhc.txt") -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Check mhc
    This function will check if the data format conforms to what our model expects

    Parameters
    ---------
    df: pd.DataFrame
        A pandas dataframe containing pairing data
    mhc_path: str
        The file path to valid mhcs 

    Returns
    ---------
    Tuple[df.DataFrame, df.DataFrame, df.DataFrame, df.DataFrame]
        Four pandas dataframe containing curated pairing data,
        pairs with peptides longer than 30,
        problematic mhca, and problematic mhcb
    """
    # input MHC that is not in the ESM dictionary and
    # whose corresponding sequence are not provided will be dropped
    with open(mhc_path, 'r') as f:
        mhc_dic_keys = f.read().splitlines()
    mhc_dic_keys = set(mhc_dic_keys)
    df_mhc_alpha_dropped = df[(~df["mhca"].isin(mhc_dic_keys)) &
                              (df['mhcaseq'] == '')].reset_index(drop=True)
    df_mhc_beta_dropped = df[(~df["mhcb"].isin(mhc_dic_keys)) &
                             (df['mhcbseq'] == '')].reset_index(drop=True)
    df = df[(df["mhca"].isin(mhc_dic_keys)) | (df['mhcaseq'] != '')]
    df = df[(df["mhcb"].isin(mhc_dic_keys)) | (df['mhcbseq'] != '')]
    df['mhca_use_seq'] = ~df['mhca'].isin(mhc_dic_keys)
    df['mhcb_use_seq'] = ~df['mhcb'].isin(mhc_dic_keys)

    df = df.reset_index(drop=True)
    print("Number of rows in processed dataset: " + str(df.shape[0]))

    return df, df_mhc_alpha_dropped, df_mhc_beta_dropped


def check_peptide(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

    # antigen peptide longer than 30 will be dropped
    df_antigen_dropped = df[df.peptide.str.len() > 30].reset_index(drop=True)
    df = df[df.peptide.str.len() <= 30].reset_index(drop=True)
    return df, df_antigen_dropped


def check_amino_acids(df_column: pd.DataFrame) -> pd.DataFrame:
    """Check amino acids are valid 
    This function checks if the amino acids in one column of a dataframe are valid amino acids

    Parameters
    ---------
    df_column: pd.DataFrame
        One column of a dataframe 

    Returns
    --------
    pd.DataFrame
        Currated column with invalid aa replaced by "_"

    """
    aa_set = set([*'ARDNCEQGHILKMFPSTWYV'])
    print('Checking amino acids...\n')
    for r in tqdm(range(df_column.shape[0])):
        wrong_aa = [aa for aa in df_column.iloc[r, 0] if aa not in aa_set]
        for aa in wrong_aa:
            df_column.iloc[r, 0] = df_column.iloc[r, 0].replace(aa, "_")
    return df_column


def check_amino_acids_columns(df: pd.DataFrame) -> pd.DataFrame:

    print("Checking if provided amino acids are valid\n")
    df['vaseq'] = check_amino_acids(df["vaseq"].to_frame())
    df['vbseq'] = check_amino_acids(df["vbseq"].to_frame())
    df['cdr3a'] = check_amino_acids(df["cdr3a"].to_frame())
    df['cdr3b'] = check_amino_acids(df["cdr3b"].to_frame())
    df['peptide'] = check_amino_acids(df["peptide"].to_frame())
    return df


def read_file(file_path: str,
              background_tcrs_dir: str = "./validation_data/",
              mhc_path: str = "./validation_data/valid_mhc.txt",
              **kwargs) -> pd.DataFrame:
    """Reads in user dataframe and performs some basic data curation 

    Parameters:
    -----------
    file_path: str
        Path to the dataframe
    background_tcrs_dir: str
        The directory with background TCR datasets 
    mhc_path: str
        The file path to valid mhcs 

    Returns
    -------
    pd.DataFrame
        A curated pandas dataframe 

    """
    print('Attempting to read in the dataframe...\n')
    df = pd.read_csv(file_path, **kwargs).fillna('')
    print("Number of rows in raw dataset: " + str(df.shape[0]))
    # Check column names
    df = check_column_names(df=df)
    # Check species
    df = check_species(df=df)
    # Check VA VB
    df, _ = check_va_vb(df=df, background_tcrs_dir=background_tcrs_dir)
    # Check MHC
    df = infer_mhc_info(df=df)
    df, _, _ = check_mhc(df=df, mhc_path=mhc_path)
    # Check peptide
    df, _ = check_peptide(df=df)
    # Check aa sequences
    df = check_amino_acids_columns(df=df)
    return df