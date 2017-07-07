import fenics
import mshr
from dolfin.cpp.function import near
from dolfin.cpp.io import File
from dolfin.cpp.mesh import Point
from fenics import *
from config_reader import Config, Geometry, Material
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


def epsilon(v: TestFunction):
    return sym(nabla_grad(v))


def sigma(u: TrialFunction, material: Material):
    E = material.young_modulus
    G = material.shear_modulus

    lambda_ = G * (E - 2 * G) / (3 * G - E)
    mu = G

    return lambda_ * nabla_div(u) * Identity(3) + 2 * mu * epsilon(u)


def solve(config: Config):
    mesh = create_mesh(config.geometry)

    vtkfile = File('results/mesh.pvd')
    vtkfile << mesh

    V = VectorFunctionSpace(mesh, 'P', 1)

    u = TrialFunction(V)
    v = TestFunction(V)

    gravity = Point(0, 0, -9.8)
    f = Constant(config.material.density * gravity)
    T = Constant((0, 0, 0))

    a = inner(sigma(u, config.material), epsilon(v)) * dx
    L = dot(f, v) * dx + dot(T, v) * ds

    bc = DirichletBC(V, Constant((0, 0, 0)),
                     lambda x, on_boundary: on_boundary and near(x[2], 0))

    u = Function(V)

    fenics.solve(a == L, u, bc)

    vtkfile = File('results/displacement.pvd')
    vtkfile << u
