from __future__ import print_function, nested_scopes, division, absolute_import

import os

from dolfin.cpp.function import near
from dolfin.cpp.io import File
from fenics import *

from mesh.mesh_factory import create_mesh


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


def output_results(u, root):
    """
    Output results in vtu format.
    :param root: base dir for results
    :param u: function to output
    """

    vtk_file = File(os.path.join(root, 'heat_conduction', '{name}.pvd'.format(name=u)))
    vtk_file << u


def main(params):
    geometry = params['geometry']

    mesh = create_mesh(geometry)

    V = FunctionSpace(mesh, 'P', 1)

    bcs = create_boundary_conditions(V, geometry['height'])

    a, L = create_variational_problem(V)

    u = Function(V, name='T')
    solve(a == L, u, bcs=bcs)

    output_results(u, params['output']['root'])
