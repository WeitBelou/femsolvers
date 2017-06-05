import os

from dolfin.cpp.io import File
from fenics import *

from boundary_conditions.factory import create_dirichlet
from config.config import Config


class StationaryHeatSolver:
    def __init__(self, config: Config, mesh: Mesh, bcs: Config):
        self._bcs = bcs
        self._config = config
        self._mesh = mesh

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
        V = FunctionSpace(self._mesh, 'P', 1)

        a, L = self.create_variational_problem(V)

        u = Function(V, name='T')
        solve(a == L, u, bcs=create_dirichlet(V, self._bcs))

        self.output_results(u, self._config['output']['root'])
