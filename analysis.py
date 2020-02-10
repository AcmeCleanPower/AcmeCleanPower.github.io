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
    adjust = kernel.sum()
    return signal.convolve2d(data, kernel, mode='same') * adjust


def generate_heatmaps_from_population(data, kernel_area=1000000, resolution=30):
    #data = np.choose(data > 0, (0, data))
    #heatmap_km = get_heatmap(data, kernel_area, resolution)
    data_points = np.choose(data > 0, (0, 1))
    heatmap_cover = get_heatmap(data_points, kernel_area, resolution)
    print heatmap_cover.max(), heatmap_cover.min()
    return heatmap_cover


def run_for_file(addr, fout1='../madagascar/points_hrsl_km2.tiff'):
    data = tiff_tools.read_array_from_tiff(addr)
    data = np.choose(data==127, (1, 0))
    print data.sum(), data.min(), data.max()
    data_points = generate_heatmaps_from_population(data, kernel_area = 1000*1000*np.pi)
    tiff_tools.save_data_derived_from_tiff(addr, data_points, fout1, dtype=np.uint16, maxcolor=8000, mincolor=0, cmap=tiff_tools.plt.cm.jet)

if __name__ == '__main__':
    import sys
    fin = sys.argv[1]
    fout1 = sys.argv[2]
    run_for_file(fin, fout1)
