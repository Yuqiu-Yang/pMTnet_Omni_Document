import os 
import sys 
from io import StringIO

import pytest 

from pMTnet_Omni_Document import __main__
import pMTnet_Omni_Document


# Version, author list, and help
def test_parser_version():
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    
    try:
        __main__.parser.parse_args(["--version"])
    except SystemExit:
        output = string_stdout.getvalue()
        expected = pMTnet_Omni_Document.__version__ + "\n"
        assert output == expected
        sys.stdout = screen_stdout
    else:
        sys.stdout = screen_stdout
        assert False

def test_parser_author():
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    
    try:
        __main__.parser.parse_args(["--author"])
    except SystemExit:
        output = string_stdout.getvalue()
        expected = pMTnet_Omni_Document.__author__ + "\n"
        assert output == expected
        sys.stdout = screen_stdout
    else:
        sys.stdout = screen_stdout
        assert False


@pytest.mark.parametrize("cmdarg", ["--help", "-h"])
def test_parser_help(cmdarg):
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    
    try:
        __main__.parser.parse_args([cmdarg])
    except SystemExit:
        output = string_stdout.getvalue()
        assert "usage: " in output
        sys.stdout = screen_stdout
    else:
        sys.stdout = screen_stdout
        assert False


@pytest.mark.parametrize("file_path, validation_data_path, output_folder_path", [
    ('./tests/test_data/test_df.csv', './validation_data', "./tests/test_data/"),
])
def test_pmtnet_omni_document_main(mocker, file_path, validation_data_path, output_folder_path):
    cmdargs: mocker.MagicMock = mocker.MagicMock() 
    cmdargs.file_path = file_path
    cmdargs.validation_data_path = validation_data_path
    cmdargs.output_folder_path = output_folder_path

    try:
        __main__.main(cmdargs=cmdargs)
    except SystemExit:
        assert True 
    except:
        assert False 


@pytest.mark.parametrize("file_path, validation_data_path, output_folder_path", [
    ('./tests/test_data/test_df.pdf', './validation_data', "./tests/test_data/"),
    ('./tests/test_data/test_df.rds', './validation_data', "./tests/test_data/"),
])
def test_pmtnet_omni_main_document_error(mocker, file_path, validation_data_path, output_folder_path):
    cmdargs: mocker.MagicMock = mocker.MagicMock() 
    cmdargs.file_path = file_path
    cmdargs.validation_data_path = validation_data_path
    cmdargs.output_folder_path = output_folder_path

    try:
        __main__.main(cmdargs=cmdargs)
    except SystemExit:
        assert False 
    except:
        assert True 

