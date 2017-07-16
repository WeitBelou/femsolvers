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

    def solve(self) -> Function:
        V = FunctionSpace(self._mesh, 'P', 1)

        a, L = self.create_variational_problem(V)

        u = Function(V, name='T')
        solve(a == L, u, bcs=create_dirichlet(V, self._bcs))
        return u
