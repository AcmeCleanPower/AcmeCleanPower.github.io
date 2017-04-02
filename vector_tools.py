import fiona
from shapely import geometry, ops


def dissolve_shapefile(fin):
    with fiona.open(fin) as fin:
        shape = map(lambda s: geometry.shape(s['geometry']), fin)
    return ops.cascaded_union(shape)


def write_dissolved_shapefile(fout, shape):
    crs = {'no_defs': True, 'ellps': 'WGS84', 'datum': 'WGS84', 'proj': 'longlat'}
    with fiona.open(fout, 'w', driver='ESRI Shapefile', schema={'geometry': 'Polygon', 'properties': {}}, crs=crs) as fout:
        fout.write({'geometry': geometry.mapping(shape), 'id': 0, 'properties': {}})
