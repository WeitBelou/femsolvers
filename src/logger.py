import logging


def configure_logger():
    """
    Configure logger and suppress some fenics logs
    """
    logging.basicConfig(format='[%(name)s] [%(levelname)s] %(asctime)s %(message)s',
                        level=logging.DEBUG)

    logging.getLogger('FFC').setLevel(logging.WARNING)
    logging.getLogger('UFL').setLevel(logging.WARNING)


def get_logger(name): return logging.getLogger(name)
