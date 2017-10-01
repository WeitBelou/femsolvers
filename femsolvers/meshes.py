import mshr
from dolfin.cpp.mesh import Mesh, Point


class UnknownGeometryType(Exception):
    def __init__(self, type_):
        self._type = type_

    def __str__(self):
        return 'Unknown geometry type: {type}'.format(type=self._type)


def create_mesh_from_geometry(geometry: dict) -> Mesh:
    """
    Creates mesh from geometry description.

    :raises: UnknownGeometryType if trying to pass geometry with unknown type
    """
    geometry_type = geometry['type']

    if geometry_type == 'cylinder':
        radius = geometry['radius']
        height = geometry['height']
        bottom = geometry['bottom']
        return _create_cylinder_mesh(radius, height, Point(bottom[0], bottom[1], bottom[2]))
    else:
        raise UnknownGeometryType(geometry_type)


def _create_cylinder_mesh(radius: float, height: float, bottom: Point = Point()) -> Mesh:
    """
    Create cylinder along Z axis.

    :param radius: radius of cylinder.
    :param height: height of cylinder, top in the positive direction of z axis.
    :param bottom: bottom point of cylinder.
    :return: mesh object
    """
    top = bottom + Point(0.0, 0.0, height)
    geometry = mshr.Cylinder(top=top, bottom=bottom, top_radius=radius, bottom_radius=radius)
    return mshr.generate_mesh(geometry, 15.0)
