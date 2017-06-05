import logging
import os
import sys

from config import parser
from meshes.factory import create_mesh
from solvers.factory import create_solver


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
        mesh = create_mesh(self._config['geometry'])

        # TODO: We need to declare function space for bcs (how to get rid of it?)
        bcs = self._config['boundary_conditions']['dirichlet']

        solver = create_solver(self._config['solver'], mesh, bcs)

        solver.run()


if __name__ == '__main__':
    runner = Runner()
    runner.run()
