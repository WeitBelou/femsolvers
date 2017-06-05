import os

from dolfin.cpp.io import File
from fenics import *

from boundary_conditions.dirichlet import create_boundary_conditions
from config.config import Config
from mesh.mesh_factory import create_mesh


class StationaryHeatSolver:
    def __init__(self, config: Config):
        self._config = config

    def create_variational_problem(self, function_space):
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

    def output_results(self, u, root):
        """
        Output results in vtu format.
        :param root: base dir for results
        :param u: function to output
        """

        vtk_file = File(os.path.join(root, 'solvers', '{name}.pvd'.format(name=u)))
        vtk_file << u

    def run(self):
        geometry = self._config['geometry']

        mesh = create_mesh(geometry)

        V = FunctionSpace(mesh, 'P', 1)

        bcs = create_boundary_conditions(V, self._config['boundary_conditions']['dirichlet'])

        a, L = self.create_variational_problem(V)

        u = Function(V, name='T')
        solve(a == L, u, bcs=bcs)

        self.output_results(u, self._config['output']['root'])
