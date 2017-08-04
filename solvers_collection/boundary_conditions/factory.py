from typing import List

from dolfin import DirichletBC, FunctionSpace

from boundary_conditions.markers import make_marker


def create_dirichlet(function_space: FunctionSpace, bcs: dict) -> List[DirichletBC]:
    """
    Create dirichlet boundary conditions
    :param function_space: function space for which boundary conditions will be created.
    :param bcs: Description of boundary conditions.
    :return: list of boundary conditions.
    """

    def _make_const_dirichlet(conf: dict):
        return DirichletBC(function_space, conf['value'], make_marker(conf['marker']))

    return list(map(_make_const_dirichlet, bcs))
