import mshr
from dolfin.cpp.mesh import Mesh, Point


def create_cylinder_mesh(radius: float, height: float, bottom: Point = Point()) -> Mesh:
    """
    Create cylinder along Z axis.

    :param radius: radius of cylinder.
    :param height: height of cylinder, top in the positive direction of z axis.
    :param bottom: bottom point of cylinder.
    :return: mesh object
    """
    top = bottom + Point(0.0, 0.0, height)
    geometry = mshr.Cylinder(top=top, bottom=bottom, top_radius=radius, bottom_radius=radius)
    return mshr.generate_mesh(geometry)
