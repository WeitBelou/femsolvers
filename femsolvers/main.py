import os
import sys
import logging

import mshr
from dolfin.cpp.mesh import Point, Mesh

import utils

BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def create_logger():
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)

    logfile = logging.FileHandler(os.path.join(BASE_DIR, 'kevlar.log'))
    logfile.setFormatter(formatter)

    logger = logging.Logger(__name__, level=logging.INFO)
    logger.addHandler(console)
    logger.addHandler(logfile)

    return logger


def create_circle(r: float) -> Mesh:
    return mshr.generate_mesh(
        mshr.Circle(Point(0.0, 0.0, 0.0), r), 100
    )


if __name__ == '__main__':
    LOG = create_logger()

    MESHES_DIR = os.path.join(BASE_DIR, 'meshes')
    RESULTS_DIR = os.path.join(BASE_DIR, 'results')

    LOG.info('Generating mesh...')
    mesh = create_circle(100.0)

    LOG.info('Writing mesh to file...')
    utils.write_data(MESHES_DIR, 'circle.pvd', mesh)
