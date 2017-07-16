from typing import List, Tuple

from fenics import *


def create_heat_variational_problem(function_space: FunctionSpace) -> Tuple:
    """
    Create linear variational problem for stationary Laplace problem.
    :param function_space:
    :return: tuple (a(u, v), L)
    """
    u = TrialFunction(function_space)
    v = TestFunction(function_space)

    a = dot(grad(u), grad(v)) * dx
    L = Constant(0.0) * v * dx

    return a, L


def solve_heat_problem(function_space: FunctionSpace, bcs: List[DirichletBC]) -> Function:
    """
    Solves heat problem on function_space with boundary conditions == bcs
    :param function_space: FunctionSpace function space on which problem has to be solved
    :param bcs: List[DirichletBC] Dirichlet boundary conditions
    :return: Function with solution of problem
    """
    a, L = create_heat_variational_problem(function_space)
    u = Function(function_space, name='T')

    solve(a == L, u, bcs=bcs)
    return u