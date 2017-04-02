import numpy as np
from scipy import signal
from skimage import morphology
import tiff_tools


def get_heatmap(data, kernel_area=1000000, resolution=30):
    kernel_radius = np.sqrt(kernel_area) / np.pi
    kernel_radius /= resolution
    if int(kernel_radius + 0.5) > int(kernel_radius):
        kernel_radius = int(kernel_radius + 0.5)
    else:
        kernel_radius = int(kernel_radius)
    kernel = morphology.disk(kernel_radius)
    adjust = np.sqrt(kernel_area) / np.pi / kernel_radius
    return signal.convolve2d(data, kernel) * adjust


def generate_heatmaps_from_population(data, kernel_area=1000000, resolution=30):
    data = np.choose(data > 0, (0, data))
    heatmap_km = get_heatmap(data, kernel_area, resolution)
    data_points = np.choose(data > 0, (0, 1))
    heatmap_cover = get_heatmap(data_points, kernel_area, resolution)
    return heatmap_km, heatmap_cover


def run_for_file(addr, fout0='../madagascar/population_hrsl_km2.tiff', fout1='../madagascar/points_hrsl_km2.tiff'):
    data = tiff_tools.read_array_from_tiff(addr)
    data_km, data_points = generate_heatmaps_from_population(data)
    tiff_tools.save_data_derived_from_tiff(addr, data_km, fout0, dtype=np.uint16, mincolor=0, cmap=tiff_tools.plt.cm.jet, logcolor=1)
    tiff_tools.save_data_derived_from_tiff(addr, data_points, fout1, dtype=np.uint16, maxcolor=255, mincolor=0, cmap=tiff_tools.plt.cm.jet)
