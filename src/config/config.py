from typing import Union

from ruamel import yaml


class Config:
    """
    Class that wraps parameters dict and give you access
    through member variables
    """

    def __init__(self, raw_params: dict):
        self._raw_params = raw_params

    def __getattribute__(self, item) -> Union[float, int, str, 'Config']:
        res = object.__getattribute__(self, '_raw_params')[item]
        if type(res) is dict:
            return Config(res)
        return res

    def __str__(self) -> str:
        return str(yaml.safe_dump(object.__getattribute__(self, '_raw_params'),
                                  default_flow_style=False))
