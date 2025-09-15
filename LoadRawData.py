"""This function loads raw data files into numpy matrices.

Input are the filenames, output the raw data matrices.
"""


import numpy as np
from aicsimageio import AICSImage


class LoadRawData:
    """Load the raw data files."""
    def __init__(self, raw_vsv_fname, raw_delta_fname):

        img_vsv         =  AICSImage(raw_vsv_fname)                                                                     # read and load file
        img_delta       =  AICSImage(raw_delta_fname)                                                                   # read and load file
        self.pix_sizeX  =  img_vsv.physical_pixel_sizes.X
        self.pix_sizeZ  =  img_vsv.physical_pixel_sizes.Z
        self.raw_vsv    =  np.rot90(np.squeeze(img_vsv.get_image_data()))[::-1]                                         # modify to have match imageJ rendering
        self.raw_delta  =  np.rot90(np.squeeze(img_delta.get_image_data()))[::-1]



