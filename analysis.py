import numpy as np
from scipy import signal
from skimage import morphology
from skimage.feature import peak_local_max
import tiff_tools
import csv


def filter_coordinates(tmp_coordinates, min_distance=33):
    coordinates = []
    cur_c = np.zeros(2)
    for c in tmp_coordinates:
        dist = np.sqrt(np.sum(np.multiply(cur_c - c, cur_c - c)))
        if dist < min_distance:
            continue
        cur_c = c
        coordinates.append(c)
    return coordinates


def get_local_maxima(data, min_distance, threshold_abs):
    tmp_coordinates = peak_local_max(data, min_distance=min_distance, threshold_abs=threshold_abs)
    coordinates = filter_coordinates(tmp_coordinates, min_distance)
    return coordinates


def write_local_maxima(fhm, fout, resolution=30, min_cutoff=100, min_distance=1000):
    distance_pixels = int(min_distance / resolution)
    data = tiff_tools.read_array_from_tiff(fhm)
    coordinates = get_local_maxima(data, distance_pixels, min_cutoff)
    latlng = tiff_tools.pixelToLatLon(fhm, [[c[1],c[0]] for c in coordinates])
    vals = [data[c[0], c[1]] for c in coordinates]
    fout = csv.writer(open(fout, 'w'))
    fout.writerow(['latitude','longitude','households'])
    for i in range(len(vals)):
        fout.writerow([latlng[i][0], latlng[i][1], vals[i]])
    return coordinates, latlng, vals


def get_heatmap(data, kernel_area=1000000, resolution=30):
    kernel_radius = np.sqrt(kernel_area) / np.pi
    kernel_radius /= resolution
    if int(kernel_radius + 0.5) > int(kernel_radius):
        kernel_radius = int(kernel_radius + 0.5)
    else:
        kernel_radius = int(kernel_radius)
    kernel = morphology.disk(kernel_radius)
    adjust = kernel.sum() * 1.0
    return signal.convolve2d(data, kernel, mode='same')


def generate_heatmaps_from_population(data, kernel_area=1000000, resolution=30):
    #data = np.choose(data > 0, (0, data))
    #heatmap_km = get_heatmap(data, kernel_area, resolution)
    data_points = np.choose(data > 0, (0, 1))
    heatmap_cover = get_heatmap(data_points, kernel_area, resolution)
    print(heatmap_cover.max(), heatmap_cover.min())
    return heatmap_cover


def run_for_file(addr, fout1='../madagascar/points_hrsl_km2.tiff'):
    data = tiff_tools.read_array_from_tiff(addr)
    data = np.choose(data==127, (1, 0))
    print(data.sum(), data.min(), data.max())
    data_points = generate_heatmaps_from_population(data, kernel_area = 1000*1000*np.pi)
    tiff_tools.save_data_derived_from_tiff(addr, data_points, fout1, dtype=np.uint16, maxcolor=8000, mincolor=0, cmap=tiff_tools.plt.cm.jet)

if __name__ == '__main__':
    import sys
    fin = sys.argv[1]
    fout1 = sys.argv[2]
    run_for_file(fin, fout1)
