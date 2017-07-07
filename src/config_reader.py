import json

from logger import get_logger

# Create logger
logger = get_logger(__name__)


class Config:
    """
    Class that encapsulates configuration.
    """

    def __init__(self, config: dict):
        self._raw_data = config

        self.material = config['material']
        self.geometry = config['geometry']
        self.drilling = config['drilling']

    def __str__(self):
        """
        Return indented data in json format.
        :return: pretty printed json representation of config.
        """
        return json.dumps(self._raw_data, indent=2)


def read_config(config_file: str) -> Config:
    with open(config_file, 'r', encoding='utf-8') as f:
        config = Config(json.load(f))
        logger.info('Config:\n%(config)s', {'config': config})
        return config
