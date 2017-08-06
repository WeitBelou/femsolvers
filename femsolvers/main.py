import os

from dolfin.cpp.io import File
from dolfin.cpp.mesh import UnitCubeMesh

from utils import fix_docker_permissions

BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

if __name__ == '__main__':
    MESHES_DIR = os.path.join(BASE_DIR, 'meshes')

    mesh_filename = os.path.join(MESHES_DIR, 'mesh.pvd')

    mesh = UnitCubeMesh(32, 32, 32)
    File(mesh_filename) << mesh

    fix_docker_permissions(MESHES_DIR)
