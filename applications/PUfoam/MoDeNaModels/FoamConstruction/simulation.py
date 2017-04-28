#!/usr/bin/env python
"""
@brief      Spatially thee-dimensional model of heat/mass transfer.
@author     Pavel Ferkl
"""
from dolfin import *
import subprocess as sp
from blessings import Terminal
def main():
    """Main function."""
    fname = "test"
    term = Terminal()
    print(
        term.yellow
        + "Working on file {}.".format(fname)
        + term.normal
    )
    # mesh_domain(fname+".geo")
    # convert_mesh(fname+".msh", fname+".xml")
    F, u, bc = preprocess(fname)
    u = integrate(F, u, bc)
    postprocess(fname, u)
    print(
        term.yellow
        + "End."
        + term.normal
    )

def mesh_domain(domain):
    """Mesh computational domain using Gmsh."""
    call = sp.Popen(['gmsh', '-3', domain])
    call.wait()

def convert_mesh(input_mesh, output_mesh):
    """Convert mesh to xml."""
    call = sp.Popen(['dolfin-convert', input_mesh, output_mesh])
    call.wait()

# class PeriodicBoundary(SubDomain):
#     """Sub domain for Periodic boundary condition."""
#     def inside(self, z, on_boundary):
#         """Top boundary is "target domain" G"""
#         return bool(z[2] < DOLFIN_EPS and z[2] > -DOLFIN_EPS and on_boundary)
#     # Map bottom boundary (H) to top boundary (G)
#     def map(self, x, y, z):
#         y[0] = x[0] - 1.0
#         y[1] = x[1]

def preprocess(fname):
    mesh = Mesh(fname+".xml")
    File(fname+"_mesh.pvd") << mesh
    boundaries = MeshFunction('size_t', mesh, fname+'_facet_region.xml')
    subdomains = MeshFunction('size_t', mesh, fname+'_physical_region.xml')
    plot(mesh)
    plot(boundaries)
    plot(subdomains)
    interactive()
    V = FunctionSpace(mesh, "CG", 2)
    # V = FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicBoundary())
    bctop = Constant(1)
    bcbot = Constant(2)
    bc = [
        DirichletBC(V, bctop, boundaries, 1),
        DirichletBC(V, bcbot, boundaries, 2)
    ]
    u = Function(V)
    v = TestFunction(V)
    D = 1
    F = inner(D*grad(u), grad(v))*dx
    return F, u, bc

def integrate(F, u, bc):
    """Integrates the equation"""
    solve(F == 0, u, bc)
    return u

def postprocess(fname, u):
    """Postprocessing of the simulation."""
    plot(u, title="Solution")
    interactive()
    file = File(fname+".pvd")
    file << u

if __name__ == "__main__":
    main()
