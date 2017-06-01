from __future__ import print_function

from mshr import Cylinder, generate_mesh

from dolfin.cpp.function import near
from dolfin.cpp.io import File
from dolfin.cpp.mesh import Point
from fenics import *


# Create cylinder mesh
def create_mesh(r, h):
    return generate_mesh(
        Cylinder(
            top=Point(0, 0, h),
            bottom=Point(0, 0, 0),
            top_radius=r, bottom_radius=r
        ), 10
    )


def create_boundary_conditions(V, h):
    u_top = Expression('0', degree=3)
    u_bottom = Expression('100', degree=3)

    bottom_bc = DirichletBC(V, u_top, lambda x, on_boundary: on_boundary and near(x[2], h))
    upper_bc = DirichletBC(V, u_bottom, lambda x, on_boundary: on_boundary and near(x[2], 0))

    return [upper_bc, bottom_bc]


if __name__ == '__main__':
    R = 10
    H = 5
    mesh = create_mesh(R, H)

    V = FunctionSpace(mesh, 'P', 1)

    # Define BC
    bcs = create_boundary_conditions(V, H)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)

    f = Constant(0.0)

    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx

    # Compute solution
    u = Function(V, name='T, K')
    solve(a == L, u, bcs=bcs)

    # Save solution in vtk file
    vtk_file = File('heat_conduction/results/solution.pvd')
    vtk_file << u
