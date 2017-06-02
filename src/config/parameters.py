class Parameters:
    """
    Class that wraps parameters dict and give you access
    through member variables
    """

    def __init__(self, raw_params: dict):
        self._raw_params = raw_params

    def __getattribute__(self, item):
        res = object.__getattribute__(self, '_raw_params')[item]
        if type(res) is dict:
            return Parameters(res)
        return res
