from __future__ import print_function

from dolfin.cpp.function import near
from dolfin.cpp.io import File
from dolfin.cpp.mesh import UnitCubeMesh
from fenics import *

# Create mesh and define function spaces
mesh = UnitCubeMesh(8, 8, 8)
V = FunctionSpace(mesh, 'P', 1)

# Define BC
u_top = Expression('0', degree=3)
u_bottom = Expression('100', degree=3)

bottom_bc = DirichletBC(V, u_top, lambda x, on_boundary: on_boundary and near(x[2], 1.0))
upper_bc = DirichletBC(V, u_bottom, lambda x, on_boundary: on_boundary and near(x[2], 0.0))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f = Constant(0.0)

a = dot(grad(u), grad(v)) * dx
L = f * v * dx

# Compute solution
u = Function(V, name='T, K')
solve(a == L, u, bcs=[upper_bc, bottom_bc])

# Save solution in vtk file
vtkfile = File('heat_conduction/results/solution.pvd')
vtkfile << u
