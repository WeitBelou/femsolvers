from ufl import Mesh

from config.config import Config
from solvers.stationary_heat_solver import solve_heat_problem


class SolverTypeNotFound(Exception):
    """
    Exception that will be raised when one try to use unknown solver
    """

    def __init__(self, solver_type: str):
        """
        Constructor from solver_type
        :param solver_type: type of solver that hasn't been founded
        """
        self.solver_type = solver_type

    def __str__(self) -> str:
        return 'Solver with type: "{type}" not found'.format(type=self.solver_type)


def create_solver(solver_type: str):
    """
    Creates solver function from solver type
    :param solver_type: string that specifies solver type
    :raise SolverTypeNotFound if there is no solver with type that was declared in config
    :return: chosen solver function
    """

    if solver_type == 'heat':
        return solve_heat_problem
    else:
        raise SolverTypeNotFound(solver_type)
