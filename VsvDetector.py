"""This function segments VSV viruses.

Input is the raw data matrix output is a labeled matrix.
"""

import numpy as np
from cellpose import models
from skimage.measure import regionprops_table
from skimage.morphology import binary_erosion, remove_small_objects


class VsvDetector:
    """Segment VSV virus."""
    def __init__(self, vsv_raw):

        mdl       =  models.CellposeModel(gpu=True, pretrained_model='nuclei')
        vsv_mask  =  mdl.eval(vsv_raw, diameter=13, flow_threshold=1.0, cellprob_threshold=0.0, stitch_threshold=0.0)[0]

        self.vsv_mask  =  vsv_mask


class VsvFilter:
    """Remove labels objects based on several parameters."""
    def __init__(self, vsv_mask_noflt, thick_thr, ratio_thr_max, ratio_thr_min, area_thr):

        vsv_mask  =  vsv_mask_noflt.copy()
        vsv_mask  =  remove_small_objects(vsv_mask, 25)
        rgp_lbls  =  regionprops_table(vsv_mask, properties=["axis_major_length", "axis_minor_length", "area", "coords"])
        aas       =  np.where(rgp_lbls["axis_minor_length"] > thick_thr)[0]
        ratios    =  rgp_lbls["axis_major_length"] / rgp_lbls["axis_minor_length"]
        bbs       =  np.where(ratios > ratio_thr_max)[0]
        ccs       =  np.where(ratios < ratio_thr_min)[0]
        dds       =  np.where(rgp_lbls["area"] > area_thr)[0]
        jjs       =  np.concatenate((aas, bbs, ccs, dds), axis=0)
        for jj in jjs:
            vsv_mask[rgp_lbls["coords"][jj][:, 0], rgp_lbls["coords"][jj][:, 1]]  =  0

        vsv_diff       =  np.sign(vsv_mask_noflt) - np.sign(vsv_mask)
        vsv_diff2show  =  vsv_diff - binary_erosion(vsv_diff)
        vsv_diff2show *=  vsv_mask_noflt

        self.vsv_mask       =  vsv_mask
        self.vsv_diff       =  vsv_diff
        self.vsv_diff2show  =  vsv_diff2show


# class VsvFilter:
#     """Filters the already segmented delta viruses."""
#     def __init__(self, vsv_mask_noflt, thickness_thr_value, up_ratio_thr_value, low_ratio_thr_value, high_area_thr_value):
#
#         vsv_mask  =  FilterThickness(vsv_mask_noflt, thickness_thr_value).img_lbls
#         vsv_mask  =  FilterRatioUpThresh(vsv_mask, up_ratio_thr_value).img_lbls
#         vsv_mask  =  FilterRatioLowThresh(vsv_mask, low_ratio_thr_value).img_lbls
#         vsv_mask  =  FilterBigAreaThreshold(vsv_mask, high_area_thr_value).img_lbls
#         vsv_diff  =  np.sign(vsv_mask_noflt) - np.sign(vsv_mask)
#         vsv_diff  =  vsv_diff - binary_erosion(vsv_diff)
#         vsv_diff *=  vsv_mask_noflt
#
#         self.vsv_mask  =  vsv_mask
#         self.vsv_diff  =  vsv_diff


# class FilterThickness:
#     """Remove labels objects based on their thickness."""
#     def __init__(self, img, thick_thr):
#
#         img_lbls  =  img.copy()
#         rgp_lbls  =  regionprops_table(img_lbls, properties=["axis_minor_length", "coords"])
#         jjs       =  np.where(rgp_lbls["axis_minor_length"] > thick_thr)[0]
#         for jj in jjs:
#             img_lbls[rgp_lbls["coords"][jj][:, 0], rgp_lbls["coords"][jj][:, 1]]  =  0
#
#         self.img_lbls  =  img_lbls
#
#
# class FilterRatioUpThresh:
#     """Remove labeled objects based on the ration between major and minor axis."""
#     def __init__(self, img, ratio_thr):
#
#         img_lbls  =  img.copy()
#         rgp_lbls  =  regionprops_table(img_lbls, properties=["axis_major_length", "axis_minor_length", "coords"])
#         ratios    =  rgp_lbls["axis_minor_length"] / rgp_lbls["axis_major_length"]
#         jjs       =  np.where(ratios > ratio_thr)[0]
#         for jj in jjs:
#             img_lbls[rgp_lbls["coords"][jj][:, 0], rgp_lbls["coords"][jj][:, 1]]  =  0
#
#         self.img_lbls  =  img_lbls


# class FilterRatioLowThresh:
#     """Remove labeled objects based on the ration between major and minor axis."""
#     def __init__(self, img, ratio_thr):
#
#         img_lbls  =  img.copy()
#         rgp_lbls  =  regionprops_table(img_lbls, properties=["axis_major_length", "axis_minor_length", "coords"])
#         ratios    =  rgp_lbls["axis_minor_length"] / rgp_lbls["axis_major_length"]
#         jjs       =  np.where(ratios < ratio_thr)[0]
#         for jj in jjs:
#             img_lbls[rgp_lbls["coords"][jj][:, 0], rgp_lbls["coords"][jj][:, 1]]  =  0
#
#         self.img_lbls  =  img_lbls


# class FilterBigAreaThreshold:
#     """Remove labeled objects with a small surface."""
#     def __init__(self, img, area_thr):
#
#         img_lbls  =  img.copy()
#         rgp_lbls  =  regionprops_table(img_lbls, properties=["area", "coords"])
#         jj2rm     =  np.where(rgp_lbls["area"] > area_thr)[0]
#         for jj in jj2rm:
#             img_lbls[rgp_lbls["coords"][jj][:, 0], rgp_lbls["coords"][jj][:, 1]]  =  0
#
#         self.img_lbls  =  img_lbls


# from skimage.color import label2rgb
# aze = np.load('/home/atrullo/Dropbox/smFisher_v5/RawData/AxialView/01072027_LM174_crop/nucs_ellips.npy')[8]
# qsd  =  np.zeros(aze.shape + (4,))
# qsd[..., :3] = label2rgb(aze)
# qsd[..., 3] = np.ones(aze.shape) * (1 - .9 * (aze==8)) * (1 - .9 * (aze==19)) * (1 - .9 * (aze==12))







