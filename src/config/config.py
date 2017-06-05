from typing import Union

from ruamel import yaml


class Config:
    """
    Class that wraps parameters dict and give you access
    through member variables
    """

    def __init__(self, raw_params: dict):
        self._raw_params = raw_params

    def __getitem__(self, item):
        return self._raw_params[item]

    def __str__(self) -> str:
        return str(yaml.safe_dump(object.__getattribute__(self, '_raw_params'),
                                  default_flow_style=False))
