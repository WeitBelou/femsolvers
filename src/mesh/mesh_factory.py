from dolfin.cpp.mesh import Point
from mshr import generate_mesh, Cylinder


class GeometryTypeNotFound(Exception):
    def __init__(self, geometry_type):
        self.message = 'Geometry with type: "{type}" not found'.format(type=geometry_type)

    def __str__(self):
        return self.message


def create_mesh(geometry):
    """
    Create cylinder mesh with bottom in (0, 0, 0).
    :type geometry: dict with geometry data
    :return: generated mesh
    """
    if geometry['type'] == 'cylinder':
        return generate_mesh(
            Cylinder(
                top=Point(0, 0, geometry['height']),
                bottom=Point(0, 0, 0),
                top_radius=geometry['radius'],
                bottom_radius=geometry['radius']
            ), geometry['resolution']
        )
    else:
        raise GeometryTypeNotFound(geometry['type'])
