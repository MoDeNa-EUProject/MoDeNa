#!/usr/bin/env python
"""
@brief      Spatially thee-dimensional model of heat/mass transfer.
@author     Pavel Ferkl
"""
from __future__ import print_function
from dolfin import *
from dolfin import near, SubDomain, File, MeshFunction, interactive
import subprocess as sp
from blessings import Terminal
XMIN = 0.0
XMAX = 1.0
YMIN = 0.0
YMAX = 1.0
ZMIN = 0.0
ZMAX = 1.0
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

class PeriodicDomain(SubDomain):
    """Class for periodic boundary conditions."""
    def inside(self, x, on_boundary):
        """
        return True if on left (XMIN) or front (YMIN) boundary AND NOT on one of
        the two slave edges
        """
        return bool((near(x[0], XMIN) or near(x[1], YMIN)) and
            (not ((near(x[0], XMAX) and near(x[1], YMIN)) or
                  (near(x[0], XMIN) and near(x[1], YMAX)))) and on_boundary)

    def map(self, x, y):
        """map y onto x"""
        if near(x[0], XMAX) and near(x[1], YMAX):
            y[0] = x[0] - (XMAX - XMIN)
            y[1] = x[1] - (YMAX - YMIN)
            y[2] = x[2]
        elif near(x[0], XMAX):
            y[0] = x[0] - (XMAX - XMIN)
            y[1] = x[1]
            y[2] = x[2]
        elif near(x[1], YMAX):
            y[0] = x[0]
            y[1] = x[1] - (YMAX - YMIN)
            y[2] = x[2]
        else:
            y[0] = -1000
            y[1] = -1000
            y[2] = -1000

def preprocess(fname):
    mesh = Mesh(fname+".xml")
    File(fname+"_mesh.pvd") << mesh
    boundaries = MeshFunction('size_t', mesh, fname+'_facet_region.xml')
    subdomains = MeshFunction('size_t', mesh, fname+'_physical_region.xml')
    plot(mesh)
    plot(boundaries)
    plot(subdomains)
    # interactive()
    V = FunctionSpace(mesh, "CG", 1)
    V = FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicDomain())
    dofmap = V.dofmap()
    print('Number of dofs', len(dofmap.dofs()))
    bctop = Constant(1)
    bcbot = Constant(2)
    bc = [
        DirichletBC(V, bctop, boundaries, 1),
        DirichletBC(V, bcbot, boundaries, 2)
    ]
    u = Function(V)
    v = TestFunction(V)
    D = Constant(1)
    # f = Expression(
    #     '1/(pow(x[0]-0.1,2)+pow(x[1]-0.1,2)+pow(x[2]-0.1,2)+1e-8)',
    #     degree=2
    # )
    F = -D*inner(grad(u), grad(v))*dx# + f*v*dx
    return F, u, bc

def integrate(F, u, bc):
    """Integrates the equation"""
    # a, L = lhs(F), rhs(F)
    # solve(a == L, u, bc)
    solve(F == 0, u, bc)
    # check periodicity
    # print(u(0.0, 0.1, 0.1))
    # print(u(1.0, 0.1, 0.1))
    # print(u(0.1, 0.0, 0.1))
    # print(u(0.1, 1.0, 0.1))
    # print(u(0.1, 0.1, 0.0))
    # print(u(0.1, 0.1, 1.0))
    return u

def postprocess(fname, u):
    """Postprocessing of the simulation."""
    plot(u, title="Solution")
    interactive()
    file = File(fname+".pvd")
    file << u

if __name__ == "__main__":
    main()
