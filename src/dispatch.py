from __future__ import print_function

import os

from heat_conduction.main import main as heat_main
from parameters import parameters_reader

if __name__ == '__main__':
    root = os.path.dirname(os.path.abspath(__file__))
    parameters_file = os.path.join(root, 'parameters.yml')

    params = parameters_reader.read_parameters(parameters_file)

    if params['solver']['type'] == 'heat':
        heat_main(params)
