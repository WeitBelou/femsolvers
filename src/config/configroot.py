from typing import List, Callable

from dolfin import FunctionSpace, DirichletBC, Function

from boundary_conditions.factory import create_dirichlet
from function_spaces import create_function_space
from meshes.factory import create_mesh
from solvers.factory import create_solver


class ConfigRoot:
    """
    Class that wraps parameters dict and give you access
    through member variables
    """

    def __init__(self, function_space: FunctionSpace, bcs: List[DirichletBC],
                 solver: Callable[[FunctionSpace, List[DirichletBC]], Function],
                 output_dir: str):
        self.function_space = function_space
        self.output_dir = output_dir
        self.solver = solver
        self.bcs = bcs

    def __str__(self) -> str:
        """
        Returns string representation of object
        :return: json-formatted string
        """
        return str({
            "function_space": self.function_space,
            "output_dir": self.output_dir,
            "solver": self.solver,
            "bcs": self.bcs,
        })


def config_root(data: dict) -> ConfigRoot:
    """
    Creates root config object from data tree
    :param data: data from config file
    :return: ConfigRoot root config for problem
    """
    mesh = create_mesh(data['geometry'])
    function_space = create_function_space(mesh, data['finite_element'])

    bcs = create_dirichlet(function_space, data['boundary_conditions']['dirichlet'])
    solver = create_solver(data['solver']['type'])

    output_dir = data['solver']['output_dir']

    return ConfigRoot(function_space, bcs, solver, output_dir)
