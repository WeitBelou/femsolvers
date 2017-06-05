import logging
import os
import sys

from dolfin.cpp.common import set_log_level

from heat_conduction.stationary_heat_solver import StationaryHeatSolver
from config import parser


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


class Runner:
    """
    Class that creates solver and passes config to it
    """

    def __init__(self):
        """
        Constructor that parses commandline args and provide default config if
        there is not other
        """
        self._configure_logger()

        parameters_file = self._get_parameters_file()

        self._config = parser.parse(parameters_file)
        self.logger.info('Config:\n%(config)s', {'config': self._config})

    def _configure_logger(self) -> None:
        """
        Configure logger and suppress some fenics logs
        """
        logging.basicConfig(format='[%(name)s] [%(levelname)s] %(asctime)s %(message)s')

        logging.getLogger('FFC').setLevel(logging.WARNING)
        logging.getLogger('UFL').setLevel(logging.WARNING)

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def _get_parameters_file(self) -> str:
        """
        Get path to config file from commandline args or return default
        :return:
        """
        if len(sys.argv) == 1:
            default_filename = 'parameters.yml'
            root = os.path.dirname(os.path.abspath(__file__))
            parameters_file = os.path.join(root, default_filename)

            self.logger.info('Using "%(params)s" config file', {'params': parameters_file})

            return parameters_file

        parameters_file = os.path.abspath(sys.argv[1])

        self.logger.info('Using "%(params)s" config file', {'params': parameters_file})

        return parameters_file

    def run(self):
        """
        Creates solver and passes config to it
        :raises SolverTypeNotFound when solver type in config invalid
        """
        solver_type = self._config['solver']['type']

        if solver_type == 'heat':
            solver = StationaryHeatSolver(config=self._config)
            solver.run()
        else:
            raise SolverTypeNotFound(solver_type)


if __name__ == '__main__':
    runner = Runner()
    runner.run()
