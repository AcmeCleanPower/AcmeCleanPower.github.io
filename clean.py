import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from skimage import morphology
from skimage.external import tifffile as tif

# borrowed from tiff_tools.py
def read_array_from_tiff(fin, band=1):
    tiff = gdal.Open(fin)
    return np.array(tiff.GetRasterBand(band).ReadAsArray())

def med_filter(mr, n=8):
	med_denoise = sp.ndimage.median_filter(mr, n)
	return med_denoise

def gauss_filter(mr, n=8):
	gauss_denoise = sp.ndimage.gaussian_filter(mr, n)
	return gauss_denoise

# def tifsave(denoised, name='denoised.tif'):
#	 tif.imsave(name, denoised)