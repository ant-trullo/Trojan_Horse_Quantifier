"""This function write in a txt file the analysis parameters.

Input are the parameters, no output, just a file written.
"""

from natsort import natsorted
import glob
import datetime
import numpy as np
from skimage.morphology import binary_dilation  # , disk
from skimage.measure import regionprops_table
from skimage.segmentation import expand_labels
import xlsxwriter
from PyQt6 import QtWidgets


def perc_overl(img_lbls1, img_lbls2):
    """Calculates the overlapping rate between two segmented images."""
    delta_on_vsv  =  np.sign(img_lbls1) * img_lbls2                                                                     # binary img_lbls1 on labeled img_lbls2 to isolate the tags
    delta_on_vsv  =  np.unique(delta_on_vsv[delta_on_vsv != 0]).size                                                    # take all the tags once and remove zeros
    return delta_on_vsv * 100 / np.unique(img_lbls1).size                                                               # ratio in percentage


class AnalysisSaver:
    """This function writes the parameters used in the analysis."""
    def __init__(self, raw_delta_fname, vsv_mask, vsv_diff, delta_mask, thickness_thr_value, up_ratio_thr_value, low_ratio_thr_value, high_area_thr_value):

        vsv4rm     =  binary_dilation(vsv_diff, footprint=np.ones((3, 3)))                                              # dilation of the vsv we remove
        rgp_delta  =  regionprops_table(delta_mask, properties=["label", "coords"])                                     # regionprops of the delta
        delta2rm   =  vsv4rm * delta_mask                                                                               # product of the labeled delta times binary dilated vsv removed
        delta2rm   =  np.unique(delta2rm[delta2rm != 0])                                                                # isolate all the tags of the overlapping without zero
        for bb in delta2rm:                                                                                             # for each delta to remove
            bb_idx  =  np.where(rgp_delta["label"] == bb)[0][0]                                                         # check the index of the tag to remove
            delta_mask[rgp_delta["coords"][bb_idx][:, 0], rgp_delta["coords"][bb_idx][:, 1]] = 0                        # remove thanks to coordinates

        if raw_delta_fname.find("RED") != -1:
            np.savez(raw_delta_fname.replace("RED", "ANALYSIS")[:-5], vsv_mask, delta_mask, np.asarray([thickness_thr_value, up_ratio_thr_value, low_ratio_thr_value, high_area_thr_value]))  # save both matrices and analmysis parameters as a list
        else:
            file_name  =  QtWidgets.QFileDialog.getSaveFileName(None, "Choose a file name and a location for it")[0]
            np.savez(file_name, vsv_mask, delta_mask, np.asarray([thickness_thr_value, up_ratio_thr_value, low_ratio_thr_value, high_area_thr_value]))  # save both matrices and analysis parameters as a list


class PostStatistics:
    """This function loads all the analysis results in a folder and writes a recap exel file."""
    def __init__(self, analyses_folder, soft_version):

        workbook  =  xlsxwriter.Workbook(analyses_folder + "/ColocRecap_Journal.xlsx", {'nan_inf_to_errors': True})

        sheet_delta  =  workbook.add_worksheet("Delta on VSV")
        sheet_vsv    =  workbook.add_worksheet("VSV on Delta")
        sheet_info   =  workbook.add_worksheet("Info")

        sheet_delta.write(0, 0, "File")
        sheet_delta.write(0, 1, "Coloc Raw")
        sheet_delta.write(0, 2, "Number of Delta")

        sheet_vsv.write(0, 0, "File")
        sheet_vsv.write(0, 1, "Coloc Raw")
        sheet_vsv.write(0, 2, "Number of VSV")

        sheet_info.write(0, 0, "Software Version")
        sheet_info.write(0, 1, soft_version)
        sheet_info.write(1, 0, "date")
        sheet_info.write(1, 1, datetime.date.today().strftime("%d%b%y"))

        files  =  natsorted(glob.glob(analyses_folder + '/*.npz'))                                                      # list of the .npz files in the folder (one per each couple of images)
        for file_cnt, file in enumerate(files):                                                                         # for each of them
            sheet_delta.write(1 + file_cnt, 0, file[file.rfind('/') + 1:-4])
            sheet_vsv.write(1 + file_cnt, 0, file[file.rfind('/') + 1:-4])
            bff_file    =  np.load(file)                                                                                # load the .npz file
            vsv_mask    =  bff_file["arr_0"]                                                                            # read the vsv mask
            sheet_vsv.write(1 + file_cnt, 2, np.unique(vsv_mask).size - 1)
            delta_mask  =  bff_file["arr_1"]                                                                            # read the delta mask
            sheet_delta.write(1 + file_cnt, 2, np.unique(delta_mask).size - 1)

            rgp_delta   =  regionprops_table(delta_mask, properties=["label", "centroid"])                              # regionprops of the delta mask
            delta_ctrs  =  np.zeros_like(delta_mask)                                                                    # initialize the matrix of the centroids of the dela
            for id_cnt, id in enumerate(rgp_delta["label"]):
                x_ctr, y_ctr              =  int(np.round(rgp_delta["centroid-0"][id_cnt])), int(np.round(rgp_delta["centroid-1"][id_cnt]))     # read and approximate as integers the centroids
                delta_ctrs[x_ctr, y_ctr]  =  id                                                                         # add a pixel with the proper tag in the centrois position

            delta_ctrs  =  expand_labels(delta_ctrs, distance=2)                                                        # expand by 2 each label (we consider that a delta is overlapping if its centroid is 2 or fewer pixels far from the delta contour)
            sheet_delta.write(1 + file_cnt, 1, perc_overl(delta_ctrs, vsv_mask))                                        # calculate the overlapping ratio of the delta centroids matrix on the vsv
            sheet_vsv.write(1 + file_cnt, 1, perc_overl(vsv_mask, delta_ctrs))                                          # as before, but vsv on delta

        workbook.close()

