from typing import Iterable

import mshr
from dolfin.cpp.mesh import Mesh, Point


class Cylinder:
    def __init__(self, radius: float, height: float, bottom: Iterable = (0.0, 0.0, 0.0)):
        """
        Create cylinder along Z axis.

        :param radius: radius of cylinder.
        :param height: height of cylinder, top in the positive direction of z axis.
        :param bottom: bottom point of cylinder.
        """
        top = Point(bottom) + Point(0.0, 0.0, height)
        geometry = mshr.Cylinder(top=top, bottom=bottom, top_radius=radius, bottom_radius=radius)

        self._mesh = mshr.generate_mesh(geometry, 15.0)

    @property
    def mesh(self) -> Mesh:
        return self._mesh
