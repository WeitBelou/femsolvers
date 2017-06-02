import os

from heat_conduction.main import main as heat_main
from parameters import parameters_reader


class SolverTypeNotFound(Exception):
    """
    Exception that will be raised when one try to use unknown solver
    """

    def __init__(self, solver_type: str):
        """
        Constructor from solver_type
        :param solver_type: type of solver that hasn't been founded
        """
        self.solver_type = solver_type

    def __str__(self) -> str:
        return 'Solver with type: "{type}" not found'.format(type=self.solver_type)


class Runner:
    """
    Class that creates solver and passes parameters to it
    """

    def __init__(self, parameters_file: str):
        """
        Constructor from path to parameters file (it has to be relative this file directory)
        :param parameters_file: relative path to file with parameters
        """
        root = os.path.dirname(os.path.abspath(__file__))
        parameters_file = os.path.join(root, parameters_file)

        self.parameters = parameters_reader.read_parameters(parameters_file)

    def run(self):
        """
        Creates solver and passes parameters to it
        :raises SolverTypeNotFound when solver type in parameters invalid
        """
        solver_type = self.parameters['solver']['type']

        if solver_type == 'heat':
            heat_main(self.parameters)
        else:
            raise SolverTypeNotFound(solver_type)


if __name__ == '__main__':
    runner = Runner('parameters.yml')
    runner.run()
