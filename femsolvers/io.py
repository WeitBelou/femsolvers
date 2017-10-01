import os
from typing import Union

from dolfin.cpp.io import File
from dolfin.cpp.mesh import Mesh
from dolfin.cpp.function import Function

BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

MESH_DIR = os.path.join(BASE_DIR, 'meshes')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')


def write_mesh(name: str, mesh: Mesh):
    """
    Writes mesh in pvd format.
    """
    _write_pvd(MESH_DIR, name, mesh)


def write_results(name: str, data: Function):
    """
    Writes results in pvd format.
    """
    _write_pvd(RESULTS_DIR, name, data)


def _write_pvd(directory: str, name: str, data: Union[Function, Mesh]):
    full_name = os.path.join(directory, name + '.pvd')

    with File(full_name) as f:
        f << data


def fix_docker_permissions(root_dir: str):
    """
    Fix docker UID issue.

    :param str root_dir:
    """
    for root, dirs, files in os.walk(root_dir):
        files_permissions = 0o666
        dirs_permissions = 0o777

        os.chmod(root, dirs_permissions)

        for directory in dirs:
            os.chmod(os.path.join(root, directory), dirs_permissions)

        for file in files:
            os.chmod(os.path.join(root, file), files_permissions)
