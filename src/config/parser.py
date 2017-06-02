import yaml

from config.parameters import Parameters


def parse(parameters_file):
    """
    Read parameter file and returns dict with config
    :param parameters_file: path to parameter file
    :return: dict with config
    """
    parameters_file = open(parameters_file, 'r')

    return Parameters(yaml.safe_load(parameters_file))
