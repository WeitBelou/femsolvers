import os
import sys

from dolfin.cpp.io import File

from config import parser
from config.configroot import ConfigRoot
from logger import get_logger, configure_logger


def get_parameters_file_path() -> str:
    """
    Get path to config file path from commandline args or return default
    :return:
    """
    logger = get_logger(__name__)

    if len(sys.argv) == 1:
        default_filename = 'config.json'
        root = os.path.dirname(os.path.abspath(__file__))
        parameters_file = os.path.join(root, default_filename)

        logger.info('Using default "%(params)s" config file', {'params': parameters_file})

        return parameters_file

    parameters_file = os.path.abspath(sys.argv[1])

    logger.info('Using "%(params)s" config file', {'params': parameters_file})

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
