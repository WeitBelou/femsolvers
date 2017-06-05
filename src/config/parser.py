from ruamel import yaml

from config.config import Config


def parse(parameters_file: str) -> Config:
    """
    Read parameter file and returns dict with config
    :param parameters_file: path to parameter file
    :return: dict with config
    """
    parameters_file = open(parameters_file, 'r')

    return Config(yaml.safe_load(parameters_file))
