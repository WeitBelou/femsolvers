from __future__ import print_function, nested_scopes, division, absolute_import

from mshr import Cylinder, generate_mesh

from dolfin.cpp.function import near
from dolfin.cpp.io import File
from dolfin.cpp.mesh import Point
from fenics import *


def read_parameters(parameters_file='parameters.yml'):
    """
    Read parameter file and returns python object with parameters
    :param parameters_file: path to parameter file
    :return:
    """
    import yaml, os

    if not os.path.isabs(parameters_file):
        root = os.path.dirname(os.path.abspath(__file__))
        parameters_file = os.path.join(root, parameters_file)

    parameters_file = open(parameters_file, 'r')
    return yaml.safe_load(parameters_file)


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
    u_top = Constant(0.0)
    u_bottom = Constant(100.0)

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
    import os
    root = os.path.dirname(os.path.abspath(__file__))

    vtk_file = File(os.path.join(root, 'results', '{name}.pvd'.format(name=u)))
    vtk_file << u


if __name__ == '__main__':
    params = read_parameters()

    geometry = params['geometry']

    mesh = create_cylinder_mesh(geometry['radius'],
                                geometry['height'])

    V = FunctionSpace(mesh, 'P', 1)

    bcs = create_boundary_conditions(V, geometry['height'])

    a, L = create_variational_problem(V)

    u = Function(V, name='T')
    solve(a == L, u, bcs=bcs)

    output_results(u)
