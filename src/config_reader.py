import json

from logger import get_logger

# Create logger
logger = get_logger(__name__)


class Config:
    """
    Class that encapsulates configuration.
    """

    def __init__(self, config):
        self.material = config.material
        self.geometry = config.geometry
        self.drilling = config.drilling


def read_config(config_file: str) -> Config:
    with open(config_file, 'r', encoding='utf-8') as f:
        config = Config(json.load(f))
        logger.info('Config:\n%(config)s', {'config': json.dumps(obj=config, indent=2)})
        return config
