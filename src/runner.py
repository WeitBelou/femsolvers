import os
import sys

from config_reader import read_config
from logger import get_logger

# Create logger
logger = get_logger(__name__)


def _get_config_file() -> str:
    """
    Get path to config file from commandline args or return default
    :return: path to config file.
    """
    parameters_file = os.path.abspath(sys.argv[1])

    logger.info('Using "%(params)s" config file', {'params': parameters_file})

    return parameters_file


def run():
    """
    Read config file and log it.
    """
    config_file = _get_config_file()
    config = read_config(config_file)


if __name__ == '__main__':
    run()
