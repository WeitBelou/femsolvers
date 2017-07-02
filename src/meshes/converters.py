from dolfin.cpp.mesh import Point


def create_point(point: list) -> Point:
    """
    Converts list of doubles that comes from json to dolfin's Point
    :param point: list from coordinates of point
    :raises: AssertionError when comes list with wrong size.
    :return: Point
    """
    assert len(point) == 3, 'List has to have 3 elements but has, {actual}'.format(actual=len(point))
    return Point(point[0], point[1], point[2])
