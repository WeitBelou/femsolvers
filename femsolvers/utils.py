import os
from typing import Union

from dolfin.cpp.io import File
from dolfin.cpp.mesh import Mesh
from dolfin.cpp.function import Function


def write_data(directory: str, name: str, data: Union[Function, Mesh]):
    """
    Output data to file ``name`` to ``directory``. Perform docker permission fix.

    :param directory: path where output file will be written.
    :param name: name of file WITH extension.
    :param data:
    """
    full_name = os.path.join(directory, name)

    with File(full_name) as f:
        f << data

    fix_docker_permissions(directory)


def fix_docker_permissions(root_dir: str):
    """
    Fix docker UID issue. NOTE: One should use this very carefully.
    This function gives permissions for EVERYONE to do EVERYTHING with all subdirs of
    root_dir.

    :param str root_dir:
    :return:
    """
    for root, dirs, files in os.walk(root_dir):
        files_permissions = 0o666
        dirs_permissions = 0o777

        os.chmod(root, dirs_permissions)

        for directory in dirs:
            os.chmod(os.path.join(root, directory), dirs_permissions)

        for file in files:
            os.chmod(os.path.join(root, file), files_permissions)
