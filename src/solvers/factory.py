from config.config import Config
from solvers.stationary_heat_solver import StationaryHeatSolver


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


def create_solver(config: Config):
    """
    Creates solver from config
    :raise SolverTypeNotFound if there is no solver with type that was declared in config
    :param config: Config with solver description
    :return: created solver
    """
    solver_type = config['type']

    if solver_type == 'heat':
        return StationaryHeatSolver(config=config)
    else:
        raise SolverTypeNotFound(solver_type)
