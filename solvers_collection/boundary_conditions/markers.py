from typing import List, Callable

from dolfin.cpp.function import near

class WrongAxis(Exception):
    def __init__(self, axis: str):
        self._axis = axis

    def __str__(self) -> str:
        return "Axis: {axis} doesn't exist".format(axis=self._axis)


class MarkerTypeNotFound(Exception):
    def __init__(self, marker_type: str):
        self._marker_type = marker_type

    def __str__(self) -> str:
        return "Marker: {marker} doesn't exist".format(marker=self._marker_type)


def make_marker(marker: dict) -> Callable[[bool, List[float]], bool]:
    marker_type = marker['type']
    if marker_type == 'plane':
        return lambda x, on_boundary: _plane(marker['center'], marker['axis'], x, on_boundary)
    raise MarkerTypeNotFound(marker_type)


def _plane(center: List[float], axis: str, x: List[float], on_boundary: bool) -> bool:
    if not on_boundary:
        return False

    axis = _get_axis_n_from_str(axis)

    return near(x[axis], center[axis])


def _get_axis_n_from_str(axis: str) -> int:
    if axis == 'x':
        return 0
    if axis == 'y':
        return 1
    if axis == 'z':
        return 2
    raise WrongAxis(axis)
