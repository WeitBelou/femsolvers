import yaml, os


def read_parameters(parameters_file):
    """
    Read parameter file and returns dict with parameters
    :param parameters_file: path to parameter file
    :return: dict with parameters
    """
    parameters_file = open(parameters_file, 'r')

    return yaml.safe_load(parameters_file)
