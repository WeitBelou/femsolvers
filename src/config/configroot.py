import json


class ConfigRoot:
    """
    Class that wraps parameters dict and give you access
    through member variables
    """

    def __init__(self, raw_params: dict):
        self._raw_params = raw_params

    def __getitem__(self, item):
        return self._raw_params[item]

    def __str__(self) -> str:
        return json.dumps(object.__getattribute__(self, '_raw_params'), indent=2)
