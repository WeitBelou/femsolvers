import logging
import os
import sys

from dolfin.cpp.common import set_log_level

from heat_conduction.main import main as heat_main
from parameters import parameters_reader


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
    Class that creates solver and passes parameters to it
    """

    def __init__(self):
        """
        Constructor that parses commandline args and provide default config if
        there is not other
        """

        # Disable debug messages from fenics
        set_log_level(logging.WARNING)

        # Create and set custom logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        logging.basicConfig(format='[%(name)s] [%(levelname)s] %(asctime)s %(message)s')

        # If there is no external parameters then use default one
        if len(sys.argv) == 1:
            default_filename = 'parameters.yml'

            root = os.path.dirname(os.path.abspath(__file__))

            self.parameters_file = os.path.join(root, default_filename)
        else:
            self.parameters_file = os.path.abspath(sys.argv[1])

        self.logger.info('Using "%(params)s" parameters file', {'params': self.parameters_file})

        self.parameters = parameters_reader.read_parameters(self.parameters_file)

    def run(self):
        """
        Creates solver and passes parameters to it
        :raises SolverTypeNotFound when solver type in parameters invalid
        """
        solver_type = self.parameters['solver']['type']

        if solver_type == 'heat':
            heat_main(self.parameters)
        else:
            raise SolverTypeNotFound(solver_type)


if __name__ == '__main__':
    runner = Runner()
    runner.run()
