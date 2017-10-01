import dolfin.cpp.mesh as mesh
import mshr


def create_circle(r: float) -> mesh.Mesh:
    return mshr.generate_mesh(
        mshr.Circle(mesh.Point(0.0, 0.0, 0.0), r), 100
    )
