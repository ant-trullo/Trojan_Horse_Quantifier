"""This function segments Delta viruses.

Input is the raw data matrix output is a labeled matrix.
"""


import numpy as np
from skimage.morphology import label, remove_small_objects
from skimage.measure import regionprops_table
from skimage.filters import gaussian, laplace, threshold_otsu
# from scipy.optimize import curve_fit


def gauss(x, mu, sigma, A):
    """Gaussian function for fitting."""
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


class DeltaDetector:
    """Segment delta viruses."""
    def __init__(self, delta_raw):

        delta_f     =  laplace(gaussian(delta_raw, 1.5))                                                  # laplacian of Gaussian filter
        thr_val     =  threshold_otsu(delta_f)
        spts3       =  delta_f > thr_val

        spts3[0, :]   =  0                                                                                              # clean borders
        spts3[-1, :]  =  0
        spts3[:, 0]   =  0
        spts3[:, -1]  =  0

        spts3       =  remove_small_objects(spts3, 4)                                                   # remove small dots, coming from segmentation
        delta_spts  =  label(spts3)                                                                               # label the remaining connected components
        rgp_delta   =  regionprops_table(delta_spts, delta_raw, properties=["mean_intensity", "coords"])
        jjs         =  np.where(rgp_delta["mean_intensity"] < 40)[0]
        for jj in jjs:
            delta_spts[rgp_delta["coords"][jj][:, 0], rgp_delta["coords"][jj][:, 1]]  =  0

        for pp in rgp_delta["coords"]:
            if np.unique(pp[:, 0]).size == 1 or np.unique(pp[:, 1]).size == 1:
                delta_spts[pp[:, 0], pp[:, 1]]  =  0

        self.delta_spts  =  delta_spts


# class DeltaDetector:
#     """Segment delta viruses."""
#     def __init__(self, delta_raw):
#
#         # delta_f     =  laplace(gaussian(delta_raw, 1.5))                                                  # laplacian of Gaussian filter
#         delta_f     =  gaussian(delta_raw, 1.2) - gaussian(delta_raw, 4.5)                              # difference of gaussian filter
#         hh          =  np.histogram(np.squeeze(delta_f.reshape((-1, 1))), bins='fd')                              # histogram of the pixel values of the filtered image
#         p_init      =  [hh[1][np.argmax(hh[0])], hh[1][46] - hh[1][44], hh[0].max()]                              # find starting points for the fitting
#         prms, cov   =  curve_fit(gauss, hh[1][:-1], hh[0], p0=p_init, maxfev=3000)                                # performa gaussian fitting on the histogram
#         spts3       =  delta_f > 10 * prms[1] + prms[0]                                                            # threshold the filtered image using the histogram fitting parameters
#         spts3       =  remove_small_objects(spts3, 4)                                                   # remove small dots, coming from segmentation
#         delta_spts  =  label(spts3)                                                                               # label the remaining connected components
#         rgp_delta   =  regionprops_table(delta_spts, delta_raw, properties=["mean_intensity", "coords"])
#         jjs         =  np.where(rgp_delta["mean_intensity"] < 40)[0]
#         for jj in jjs:
#             delta_spts[rgp_delta["coords"][jj][:, 0], rgp_delta["coords"][jj][:, 1]]  =  0
#
#         self.delta_spts        =  delta_spts

