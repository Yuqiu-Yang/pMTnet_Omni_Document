"""
pMTnet Omni Document

Document for pMTnet Omni
"""

import os 

__version__ = "0.0.5"
__author__ = "Yi Han, Yuqiu Yang, Tao Wang"

# location would be pMTnet_Omni_Document/pMTnet_Omni_Document
location = os.path.dirname(os.path.realpath(__file__))
# We join .. to go back a level where the validation_data folder resides
if "pMTnet_Omni_Document" in os.path.dirname(location):
    validation_data_path = os.path.join(location, "..", "validation_data")
else:
    validation_data_path = os.path.join(location, "validation_data")
