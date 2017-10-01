import os

from femsolvers import utils

BASE_DIR = os.path.dirname(os.path.realpath(__file__))
utils.fix_docker_permissions(os.path.join(BASE_DIR, '.cache'))
