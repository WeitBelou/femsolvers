# Configure logger
import logging


def setup_logger() -> None:
    logging.basicConfig(format='[%(name)s] [%(levelname)s] %(asctime)s %(message)s',
                        level=logging.INFO)
    logging.getLogger('FFC').setLevel(logging.WARNING)
    logging.getLogger('UFL').setLevel(logging.WARNING)
    pass


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    return logger
