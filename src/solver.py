import mshr
from dolfin.cpp.io import File
from dolfin.cpp.mesh import Point
from fenics import Mesh

from config_reader import Config, Geometry
from logger import get_logger


def create_mesh(geometry: Geometry) -> Mesh:
    radius = geometry.diameter / 2
    full_height = geometry.full_height
    mesh = mshr.generate_mesh(
        mshr.Cylinder(top=Point(0.0, 0.0, full_height),
                      bottom=Point(0.0, 0.0, 0.0),
                      top_radius=radius,
                      bottom_radius=radius),
        128  # Resolution
    )
    get_logger(__name__).info('Created mesh: %(mesh)s', {'mesh': mesh})
    return mesh


def solve(config: Config):
    mesh = create_mesh(config.geometry)

    vtkfile = File('results/mesh.pvd')
    vtkfile << mesh
