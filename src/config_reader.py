import json

from dolfin.cpp.mesh import Point

from logger import get_logger

# Create logger
logger = get_logger(__name__)


class BaseConfig:
    """
    Base class for all configuration classes
    """

    def __init__(self, config: dict):
        """
        Assigns _raw_data to use in all derived classes
        :param config: dict with config data
        """
        self._raw_data = config

    def __str__(self):
        """
        Return indented data in json format.
        :return: pretty printed json representation of config.
        """
        return json.dumps(self._raw_data, indent=2)


class Material(BaseConfig):
    """
    Simple wrapper for material config.
    """

    def __init__(self, config: dict):
        super(Material, self).__init__(config)

    @property
    def density(self) -> float:
        return self._raw_data['density']

    @property
    def young_modulus(self) -> float:
        return self._raw_data['young_modulus']

    @property
    def shear_modulus(self) -> float:
        return self._raw_data['shear_modulus']

    @property
    def thermal_diffusivity(self) -> float:
        return self._raw_data['thermal_diffusivity']

    @property
    def thermal_expansion(self) -> float:
        return self._raw_data['thermal_expansion']


class Geometry(BaseConfig):
    """
    Simple wrapper for geometry config.
    """

    def __init__(self, config: dict):
        super(Geometry, self).__init__(config)

    @property
    def diameter(self) -> float:
        return self._raw_data['diameter']

    @property
    def full_height(self) -> float:
        return self._raw_data['full_height']

    @property
    def above_water_height(self) -> float:
        return self._raw_data['above_water_height']


class Drilling(BaseConfig):
    """
    Simple wrapper for drilling config.
    """

    def __init__(self, config: dict):
        super(Drilling, self).__init__(config)

    @property
    def shape(self) -> str:
        return self._raw_data['shape']

    @property
    def diameter(self) -> float:
        return self._raw_data['diameter']

    @property
    def center(self) -> Point:
        return Point(self._raw_data['center'])

    @property
    def mass(self) -> float:
        return self._raw_data['mass']


class Config(BaseConfig):
    """
    Simple wrapper for root config.
    """

    def __init__(self, config: dict):
        super(Config, self).__init__(config)

    @property
    def material(self):
        return Material(self._raw_data['material'])

    @property
    def geometry(self):
        return Geometry(self._raw_data['geometry'])

    @property
    def drilling(self):
        return Drilling(self._raw_data['drilling'])


def read_config(config_file: str) -> Config:
    with open(config_file, 'r', encoding='utf-8') as f:
        config = Config(json.load(f))
        logger.info('Config:\n%(config)s', {'config': config})
        return config
