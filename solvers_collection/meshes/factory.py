from dolfin.cpp.mesh import Point
from mshr import generate_mesh, Cylinder, Box

from ..logger import get_logger
from .converters import create_point


class GeometryTypeNotFound(Exception):
    def __init__(self, geometry_type: str):
        self.geometry_type = geometry_type

    def __str__(self):
        return 'Geometry with type: "{type}" not found'.format(type=self.geometry_type)


def create_mesh(geometry: dict):
    """
    Create mesh from geometry data.
    :type geometry: dict with geometry data
    :return: generated mesh.
    """
    if geometry['type'] == 'cylinder':
        get_logger(__name__).info('Using mesh with type: {geometry_type}'.format(
            geometry_type=geometry['type']
        ))

        return generate_mesh(
            Cylinder(
                top=Point(0, 0, geometry['height']),
                bottom=Point(0, 0, 0),
                top_radius=geometry['radius'],
                bottom_radius=geometry['radius']
            ), geometry['resolution']
        )
    if geometry['type'] == 'box':
        return generate_mesh(
            Box(a=Point(create_point(geometry['a'])),
                b=Point(create_point(geometry['b']))
                ), geometry['resolution']
        )
    else:
        raise GeometryTypeNotFound(geometry['type'])
