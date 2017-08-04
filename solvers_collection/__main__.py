import os
import sys

from dolfin.cpp.io import File

from solvers_collection.config import parser
from solvers_collection.config.configroot import ConfigRoot
from solvers_collection.logger import get_logger, configure_logger

LOGGER = get_logger(__name__)


def get_parameters_file_path() -> str:
    """
    Get path to config file path from commandline args or return default
    :return:
    """

    parameters_file = os.path.abspath(sys.argv[1])

    LOGGER.info('Using "{}" config file'.format(parameters_file))

    return parameters_file


def solve(config: ConfigRoot):
    """
    Creates solver and passes config to it
    :raises SolverTypeNotFound when solver type in config invalid
    """
    solution = config.solver(config.function_space, config.bcs)

    vtkfile = File(os.path.join(config.output_dir, 'result.pvd'))
    vtkfile << solution


if __name__ == '__main__':
    configure_logger()

    parameters_file = get_parameters_file_path()

    config = parser.parse(parameters_file)

    solve(config)
