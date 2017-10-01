import os
import stat

from femsolvers import utils

BASE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'util_test_data')

TEST_DIRS = [
    os.path.join(BASE_DIR, 'meshes'),
]

TEST_FILES = [
    os.path.join(BASE_DIR, 'meshes', 'mesh.msh')
]


def setup_module():
    os.makedirs(*TEST_DIRS)
    for file in TEST_FILES:
        with open(file, mode='x') as f:
            pass


def teardown_module():
    for file in TEST_FILES:
        os.remove(file)

    os.removedirs(*TEST_DIRS)


def test_docker_fix():
    utils.fix_docker_permissions(BASE_DIR)

    for directory in TEST_DIRS:
        assert stat.S_IMODE(os.stat(directory).st_mode) == 0o777

    for file in TEST_FILES:
        assert stat.S_IMODE(os.stat(file).st_mode) == 0o666
