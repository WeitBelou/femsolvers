import os

from dolfin.cpp.mesh import UnitCubeMesh

from utils import write_data

BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

if __name__ == '__main__':
    MESHES_DIR = os.path.join(BASE_DIR, 'meshes')

    mesh = UnitCubeMesh(32, 32, 32)
    write_data(MESHES_DIR, 'mesh.pvd', mesh)
