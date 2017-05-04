#!/usr/bin/env python
"""
@brief      Spatially thee-dimensional model of heat/mass transfer.
@author     Pavel Ferkl
"""
from __future__ import print_function
import subprocess as sp
import fenics as fe
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
    F, u, bcond = preprocess(fname)
    sol = integrate(F, u, bcond)
    postprocess(fname, sol)
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

class PeriodicDomain(fe.SubDomain):
    """Class for periodic boundary conditions."""
    def inside(self, x, on_boundary):
        """
        return True if on left (XMIN) or front (YMIN) boundary AND NOT on one of
        the two slave edges
        """
        res = bool(
            (fe.near(x[0], XMIN) or fe.near(x[1], YMIN)) and not (
                (fe.near(x[0], XMAX) and fe.near(x[1], YMIN)) or
                (fe.near(x[0], XMIN) and fe.near(x[1], YMAX))
            ) and on_boundary)
        return res

    def map(self, x, y):
        """map y onto x"""
        if fe.near(x[0], XMAX) and fe.near(x[1], YMAX):
            y[0] = x[0] - (XMAX - XMIN)
            y[1] = x[1] - (YMAX - YMIN)
            y[2] = x[2]
        elif fe.near(x[0], XMAX):
            y[0] = x[0] - (XMAX - XMIN)
            y[1] = x[1]
            y[2] = x[2]
        elif fe.near(x[1], YMAX):
            y[0] = x[0]
            y[1] = x[1] - (YMAX - YMIN)
            y[2] = x[2]
        else:
            y[0] = -1000
            y[1] = -1000
            y[2] = -1000

def preprocess(fname):
    mesh = fe.Mesh(fname+".xml")
    fe.File(fname+"_mesh.pvd") << mesh
    boundaries = fe.MeshFunction('size_t', mesh, fname+'_facet_region.xml')
    subdomains = fe.MeshFunction('size_t', mesh, fname+'_physical_region.xml')
    fe.plot(mesh)
    fe.plot(boundaries)
    fe.plot(subdomains)
    # interactive()
    # V = FunctionSpace(mesh, "CG", 1)
    V = fe.FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicDomain())
    dofmap = V.dofmap()
    print('Number of dofs', len(dofmap.dofs()))
    bctop = fe.Constant(1)
    bcbot = fe.Constant(2)
    bc = [
        fe.DirichletBC(V, bctop, boundaries, 1),
        fe.DirichletBC(V, bcbot, boundaries, 2)
    ]
    u = fe.Function(V)
    v = fe.TestFunction(V)
    D = fe.Constant(1)
    f = fe.Expression(
        '1/(pow(x[0]-0.1,2)+pow(x[1]-0.1,2)+pow(x[2]-0.1,2)+1e-8)',
        degree=2
    )
    F = -D*fe.inner(fe.grad(u), fe.grad(v))*fe.dx + f*v*fe.dx
    return F, u, bc

def integrate(F, u, bc):
    """Integrates the equation"""
    # a, L = lhs(F), rhs(F)
    # solve(a == L, u, bc)
    fe.solve(F == 0, u, bc)
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
    fe.plot(u, title="Solution")
    fe.interactive()
    fe.File(fname+".pvd") << u

if __name__ == "__main__":
    main()
