"""
pMTnet Omni Document

Document for pMTnet Omni
"""

import os 

__version__ = "0.0.15"
__author__ = "Yi Han, Yuqiu Yang, Tao Wang"

# location would be pMTnet_Omni_Document/pMTnet_Omni_Document
location = os.path.dirname(os.path.realpath(__file__))

validation_data_path = os.path.join(location, "validation_data")

from pMTnet_Omni_Document.data_curation import read_file
