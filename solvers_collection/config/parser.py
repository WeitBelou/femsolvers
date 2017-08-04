import json
from json import JSONDecodeError

from config.configroot import ConfigRoot, config_root
from logger import get_logger


def parse(parameters_file: str) -> ConfigRoot:
    """
    Read parameter file and returns dict with config
    :param parameters_file: path to parameter file
    :return: dict with config
    """
    with open(parameters_file, 'r') as parameters:
        try:
            config = json.load(parameters)
            get_logger(__name__).debug(
                'Loaded config:\n%(config)s',
                {'config': json.dumps(config, indent=1)})

            config = config_root(config)
            get_logger(__name__).debug('Parsed config:\n%(config)s', {'config': config})
            return config
        except JSONDecodeError as e:
            get_logger(__name__).error(
                'Can not load json: {parameters_file}.\nReason: {msg}'.format(
                    parameters_file=parameters_file, msg=e.msg
                ))
            exit(1)
