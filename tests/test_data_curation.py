import pytest
import pandas as pd
from copy import deepcopy

from pMTnet_Omni_Document.data_curation import check_column_names,\
    check_species,\
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
########################
# Column names


@pytest.mark.parametrize("df, expected, success_assertion", [
    (df_0, common_column_names, True),
    (df_0.drop(columns=["va", "vb", 'mhc']), common_column_names, True),
    (df_0.drop(columns=["vaseq", "vbseq", 'mhcseq']),
     common_column_names, True),
    (df_0.drop(columns=['cdr3a']), None, False),
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


@pytest.mark.parametrize("df, expected", [
    (df_0, 0)
])
def test_check_va_vb(df, expected):
    df, invalid_df = check_va_vb(df, background_tcrs_dir="./validation_data/")
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
        df, _, _ = check_mhc(df, mhc_path="./validation_data/valid_mhc.txt")
        assert True
    except:
        assert False

########################
# encode mhc sequences 


@pytest.mark.parametrize("df, output_path", [
    (df_0, "./tests/test_data/test_df_mhc_seq_dict.pickle")
])
def test_encode_mhc_seq(df, output_path):
    try:
        df = check_column_names(df)
        df = infer_mhc_info(df)
        df, _, _ = check_mhc(df, mhc_path="./validation_data/valid_mhc.txt")
        df['mhca_use_seq'] = [True for i in range(4)]
        encode_mhc_seq(df=df, output_path=output_path)
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
# read_file
@pytest.mark.parametrize("", [
    ()
])
def test_read_file():
    try:
        df = read_file(file_path='./tests/test_data/test_df.csv',
                       background_tcrs_dir="./validation_data/",
                       mhc_path="./validation_data/valid_mhc.txt",
                       output_path=None,
                       sep=",",
                       header=0)
        assert True
    except:
        assert False
