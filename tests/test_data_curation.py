import pytest
import os
import pandas as pd
from copy import deepcopy

from pMTnet_Omni_Document.data_curation import check_column_names,\
    check_species,\
    check_v_gene_allele,\
    check_va_vb,\
    infer_mhc_info,\
    check_mhc,\
    encode_mhc_seq,\
    check_peptide,\
    check_amino_acids_columns,\
    read_file

df_0 = pd.read_csv('./tests/test_data/test_df.csv', sep=',')
common_column_names = ["va", "cdr3a", "vaseq", "vb", "cdr3b", "vbseq",
                       "peptide", "mhc", "mhcseq",
                       "tcr_species", "pmhc_species"]

background_tcrs_dir = "./pMTnet_Omni_Document/validation_data"
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

########################
# Column names


@pytest.mark.parametrize("df, expected, success_assertion", [
    (df_0, common_column_names, True),
    (df_0.drop(columns=["va", "vb", 'mhc']), common_column_names, True),
    (df_0.drop(columns=["vaseq", "vbseq", 'mhcseq']),
     common_column_names, True),
    (df_0.drop(columns=['cdr3a']), None, True),
    (df_0.drop(columns=['mhc', 'mhcseq']), None, False),
    (df_0.drop(columns=['va', 'vaseq']), None, False),
    (df_0.drop(columns=['vb', 'vbseq']), None, False),
])
def test_check_column_names(df, expected, success_assertion):
    try:
        df = check_column_names(df)
        if expected is not None:
            assert set(df.columns) == set(expected)
        else:
            assert success_assertion
    except:
        assert not success_assertion


########################
# Species
@pytest.mark.parametrize("tcr_species, pmhc_species", [
    (['ape', 'ape', 'rat', 'rat'], None),
    (None, ['ape', 'ape', 'rat', 'rat'])
])
def test_check_species(tcr_species, pmhc_species):
    df = deepcopy(df_0)
    if tcr_species is not None:
        df['tcr_species'] = tcr_species
    if pmhc_species is not None:
        df['pmhc_species'] = pmhc_species
    try:
        check_species(df)
        assert False
    except:
        assert True

########################
# va vb


@pytest.mark.parametrize("a_reference_df, b_reference_df", [
    (human_alpha_tcrs, human_beta_tcrs)
])
def test_check_v_gene_allele(a_reference_df, b_reference_df):
    df = deepcopy(df_0)
    df.at[0, 'va'] = 'TRAV19'
    try:
        check_v_gene_allele(df, a_reference_df, b_reference_df)
        assert True
    except:
        assert False


@pytest.mark.parametrize("df, expected", [
    (df_0, 0)
])
def test_check_va_vb(df, expected):
    df, invalid_df = check_va_vb(df, background_tcrs_dir="./pMTnet_Omni_Document/validation_data/")
    assert invalid_df.shape[0] == expected

########################
# mhc_info


@pytest.mark.parametrize("mhc, mhcseq, expected, success_assertion", [
    (['A', 'DR', 'D', 'IA'], ['AA', 'AA', 'AA', 'AA/AA'], ['human class i',
     'human class ii', 'mouse class i', 'mouse class ii'], True),
    (['DR/DR', '', 'IA', ''], ['AA/AA', 'AA/AA', 'AA/AA', 'AA/AA'],
     ['human class ii', 'human', 'mouse class ii', 'mouse'], True),
    (['', 'Z', 'D', 'IA'], ['', '', 'AA', 'AA/AA'], None, False),
    (['A', 'DR', '', 'Z'], ['AA', 'AA', '', ''], None, False),
])
def test_infer_mhc_info(mhc, mhcseq, expected, success_assertion):
    df = deepcopy(df_0)
    if mhc is not None:
        df['mhc'] = mhc
    if mhcseq is not None:
        df['mhcseq'] = mhcseq
    try:
        df = infer_mhc_info(df)
        if expected is not None:
            assert all(df['mhc_class'] == expected)
        else:
            assert success_assertion
    except:
        assert not success_assertion

########################
# mhc


@pytest.mark.parametrize("df", [
    (df_0)
])
def test_check_mhc(df):
    try:
        df = check_column_names(df)
        df = infer_mhc_info(df)
        df, _, _ = check_mhc(df, mhc_path="./pMTnet_Omni_Document/validation_data/valid_mhc.txt")
        assert True
    except:
        assert False

########################
# peptide


@pytest.mark.parametrize("df", [
    (df_0)
])
def test_check_peptide(df):
    try:
        df, _ = check_peptide(df)
        assert True
    except:
        assert False

########################
# aa columns


@pytest.mark.parametrize("df", [
    (df_0)
])
def test_check_amino_acids_columns(df):
    df = check_amino_acids_columns(df)
    assert True


########################
# encode mhc sequences


@pytest.mark.parametrize("df", [
    (df_0)
])
def test_encode_mhc_seq(df):
    try:
        df = check_column_names(df)
        df = infer_mhc_info(df)
        df, _, _ = check_mhc(df, mhc_path="./pMTnet_Omni_Document/validation_data/valid_mhc.txt")
        df['mhcaseq'] = ['AA' for i in range(4)]
        df['mhca_use_seq'] = [True for i in range(4)]
        encode_mhc_seq(df=df)
        assert True
    except:
        assert False


########################
# read_file
@pytest.mark.parametrize("", [
    ()
])
def test_read_file():
    try:
        df, mhc_seq_dict = read_file(file_path='./tests/test_data/test_df.csv',
                                     save_results=True,
                                     output_folder_path=None,
                                     sep=",",
                                     header=0)
        try:
            os.remove("./tests/test_data/df_curated.csv")
            os.remove("./tests/test_data/df_curated_antigen_dropped.csv")
            os.remove("./tests/test_data/df_curated_invalid_v.csv")
            os.remove("./tests/test_data/df_curated_mhc_alpha_dropped.csv")
            os.remove("./tests/test_data/df_curated_mhc_beta_dropped.csv")
            os.remove("./tests/test_data/mhc_seq_dict.json")
        except:
            pass
        assert True
    except:
        assert False
