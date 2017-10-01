import logging
import os

from femsolvers.io import BASE_DIR


def create():
    formatter = logging.Formatter('[%(levelname)s] %(asctime)s "%(message)s"')

    console = logging.StreamHandler()
    console.setFormatter(formatter)

    logfile = logging.FileHandler(os.path.join(BASE_DIR, 'kevlar.log'))
    logfile.setFormatter(formatter)

    logger = logging.Logger(__name__, level=logging.INFO)
    logger.addHandler(console)
    logger.addHandler(logfile)

    return logger
