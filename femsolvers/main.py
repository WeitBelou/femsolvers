import os

from meshes import create_mesh_from_geometry
from utils import write_data

BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

if __name__ == '__main__':
    MESHES_DIR = os.path.join(BASE_DIR, 'meshes')

    mesh = create_mesh_from_geometry(dict(
        type='cylinder', radius=150.0, height=10.0, bottom=[0.0, 0.0, 0.0]
    ))
    write_data(MESHES_DIR, 'mesh.pvd', mesh)
