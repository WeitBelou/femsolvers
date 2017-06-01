from __future__ import print_function

from mshr import Cylinder, generate_mesh

from dolfin.cpp.function import near
from dolfin.cpp.io import File
from dolfin.cpp.mesh import Point
from fenics import *


def create_cylinder_mesh(r, h):
    """
    Create cylinder mesh with bottom in (0, 0, 0).
    :param r: cylinder radius
    :param h: cylinder height
    :return: generated mesh
    """
    return generate_mesh(
        Cylinder(
            top=Point(0, 0, h),
            bottom=Point(0, 0, 0),
            top_radius=r, bottom_radius=r
        ), 10
    )


def create_boundary_conditions(function_space, h):
    """
    Create boundary conditions: two planes with fixed temperature.
    :param function_space: function space for which boundary conditions will be created.
    :param h: distance between planes with fixed temperature.
    :return: list of boundary conditions.
    """
    u_top = Expression('0', degree=3)
    u_bottom = Expression('100', degree=3)

    bottom_bc = DirichletBC(function_space, u_top, lambda x, on_boundary: on_boundary and near(x[2], h))
    upper_bc = DirichletBC(function_space, u_bottom, lambda x, on_boundary: on_boundary and near(x[2], 0))

    return [upper_bc, bottom_bc]


def create_variational_problem(function_space):
    """
    Create linear variational problem for stationary Laplace problem.
    :param function_space:
    :return: pair (a(u, v), L)
    """
    u = TrialFunction(function_space)
    v = TestFunction(function_space)

    a = dot(grad(u), grad(v)) * dx
    L = Constant(0.0) * v * dx

    return a, L


def output_results(u):
    """
    Output results in vtu format.
    :param u: function to output
    """
    vtk_file = File('heat_conduction/results/{name}.pvd'.format(name=u))
    vtk_file << u


if __name__ == '__main__':
    R = 10
    H = 5

    mesh = create_cylinder_mesh(R, H)

    V = FunctionSpace(mesh, 'P', 1)

    bcs = create_boundary_conditions(V, H)

    a, L = create_variational_problem(V)

    u = Function(V, name='T')
    solve(a == L, u, bcs=bcs)

    output_results(u)
