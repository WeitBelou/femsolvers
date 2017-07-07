# Configure logger
import logging


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(format='[%(name)s] [%(levelname)s] %(asctime)s %(message)s',
                        level=logging.INFO)

    logging.getLogger('FFC').setLevel(logging.WARNING)
    logging.getLogger('UFL').setLevel(logging.WARNING)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    return logger
