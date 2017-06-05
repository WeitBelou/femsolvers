import ruamel.yaml as yaml

from config.config import Config


def parse(parameters_file):
    """
    Read parameter file and returns dict with config
    :param parameters_file: path to parameter file
    :return: dict with config
    """
    parameters_file = open(parameters_file, 'r')

    return Config(yaml.safe_load(parameters_file))
