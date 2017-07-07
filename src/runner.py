import json
import os
import sys

from logger import get_logger

# Create logger
logger = get_logger(__name__)


def _get_parameters_file() -> str:
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
    parameters_file = _get_parameters_file()

    with open(parameters_file, 'r', encoding='utf-8') as f:
        config = json.load(f)
        logger.info('Config:\n%(config)s', {'config': json.dumps(obj=config, indent=2)})


if __name__ == '__main__':
    run()
