import json
from json import JSONDecodeError

from config.configroot import ConfigRoot
from logger import get_logger


def parse(parameters_file: str) -> ConfigRoot:
    """
    Read parameter file and returns dict with config
    :param parameters_file: path to parameter file
    :return: dict with config
    """
    with open(parameters_file, 'r') as parameters:
        try:
            config = ConfigRoot(json.load(parameters))
            get_logger(__name__).debug('Config parsed:\n%(config)s', {'config': config})
            return config
        except JSONDecodeError as e:
            get_logger(__name__).error(
                'Can not load json: {parameters_file}.\nReason: {msg}'.format(
                    parameters_file=parameters_file, msg=e.msg
                ))
            exit(1)
