#!/usr/bin/env python
"""Spatially three-dimensional model of heat/mass transfer.

Usage:
    simulation.py [-h | --help] [-i input_file] [--verbose]

Options:
    -h --help       Show this screen.
    -i input_file   Json file with inputs.
    --verbose       Print more information.

@author: Pavel Ferkl
"""
from __future__ import print_function
import json
import subprocess as sp
import fenics as fe
from blessings import Terminal
from docopt import docopt
XMIN = 0.0
XMAX = 1.0
YMIN = 0.0
YMAX = 1.0
ZMIN = 0.0
ZMAX = 1.0
def main():
    """Main function. Organizes workflow."""
    fname = str(INPUTS['filename'])
    term = Terminal()
    print(
        term.yellow
        + "Working on file {}.".format(fname)
        + term.normal
    )
    if INPUTS['mesh_domain']:
        mesh_domain(fname+".geo")
    if INPUTS['convert_mesh']:
        convert_mesh(fname+".msh", fname+".xml")
    system_matrix, field, bcs = preprocess(fname)
    field = integrate(system_matrix, field, bcs)
    postprocess(fname, field)
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
    """Convert mesh to xml using dolfin-convert."""
    call = sp.Popen(['dolfin-convert', input_mesh, output_mesh])
    call.wait()

class PeriodicDomain(fe.SubDomain):
    """Class for periodic boundary conditions."""
    def inside(self, pos, on_boundary):
        """
        return True if on left (XMIN) or front (YMIN) boundary AND NOT on one of
        the two slave edges
        """
        res = bool(
            (fe.near(pos[0], XMIN) or fe.near(pos[1], YMIN)) and not (
                (fe.near(pos[0], XMAX) and fe.near(pos[1], YMIN)) or
                (fe.near(pos[0], XMIN) and fe.near(pos[1], YMAX))
            ) and on_boundary)
        return res

    def map(self, pos_a, pos_b):
        """map pos_b onto pos_a"""
        if fe.near(pos_a[0], XMAX) and fe.near(pos_a[1], YMAX):
            pos_b[0] = pos_a[0] - (XMAX - XMIN)
            pos_b[1] = pos_a[1] - (YMAX - YMIN)
            pos_b[2] = pos_a[2]
        elif fe.near(pos_a[0], XMAX):
            pos_b[0] = pos_a[0] - (XMAX - XMIN)
            pos_b[1] = pos_a[1]
            pos_b[2] = pos_a[2]
        elif fe.near(pos_a[1], YMAX):
            pos_b[0] = pos_a[0]
            pos_b[1] = pos_a[1] - (YMAX - YMIN)
            pos_b[2] = pos_a[2]
        else:
            pos_b[0] = -1000
            pos_b[1] = -1000
            pos_b[2] = -1000

def preprocess(fname):
    """Loads mesh, defines system of equations and prepares system matrix."""
    mesh = fe.Mesh(fname+".xml")
    if INPUTS['saving']['mesh']:
        fe.File(fname+"_mesh.pvd") << mesh
    boundaries = fe.MeshFunction('size_t', mesh, fname+'_facet_region.xml')
    if INPUTS['saving']['boundaries']:
        fe.File(fname+"_subdomains.pvd") << boundaries
    subdomains = fe.MeshFunction('size_t', mesh, fname+'_physical_region.xml')
    if INPUTS['saving']['subdomains']:
        fe.File(fname+"_subdomains.pvd") << subdomains
    if INPUTS['plotting']['mesh']:
        fe.plot(mesh, title='Mesh')
    if INPUTS['plotting']['boundaries']:
        fe.plot(boundaries, title='Boundaries')
    if INPUTS['plotting']['subdomains']:
        fe.plot(subdomains, title='Subdomains')
    fun_space = fe.FunctionSpace(
        mesh,
        "CG",
        INPUTS['element_degree'],
        constrained_domain=PeriodicDomain()
    )
    if ARGS['--verbose']:
        dofmap = fun_space.dofmap()
        print('Number of DOFs:', len(dofmap.dofs()))
    field = fe.TrialFunction(fun_space)
    test_func = fe.TestFunction(fun_space)
    dx = fe.Measure('dx', domain=mesh, subdomain_data=subdomains)
    cond1 = fe.Constant(INPUTS['conductivity']['gas'])
    cond2 = fe.Constant(INPUTS['conductivity']['solid'])
    system_matrix = (
        -cond1*fe.inner(fe.grad(field), fe.grad(test_func))*dx(1)
        -cond2*fe.inner(fe.grad(field), fe.grad(test_func))*dx(2)
    )
    bctop = fe.Constant(INPUTS['boundary_conditions']['top'])
    bcbot = fe.Constant(INPUTS['boundary_conditions']['bottom'])
    bcs = [
        fe.DirichletBC(fun_space, bctop, boundaries, 1),
        fe.DirichletBC(fun_space, bcbot, boundaries, 2)
    ]
    field = fe.Function(fun_space)
    return system_matrix, field, bcs

def integrate(system_matrix, field, bcs):
    """Integrates the equation"""
    left_side, right_side = fe.lhs(system_matrix), fe.rhs(system_matrix)
    fe.solve(left_side == right_side, field, bcs)
    if ARGS['--verbose']:
        print('Checking periodicity:')
        print('Value at XMIN:', field(XMIN, (YMIN + YMAX)/3, (ZMIN + ZMAX)/3))
        print('Value at XMAX:', field(XMAX, (YMIN + YMAX)/3, (ZMIN + ZMAX)/3))
        print('Value at YMIN:', field((XMIN + XMAX)/3, YMIN, (ZMIN + ZMAX)/3))
        print('Value at YMAX:', field((XMIN + XMAX)/3, YMAX, (ZMIN + ZMAX)/3))
        print('Value at ZMIN:', field((XMIN + XMAX)/3, (YMIN + YMAX)/3, ZMIN))
        print('Value at ZMAX:', field((XMIN + XMAX)/3, (YMIN + YMAX)/3, ZMAX))
    return field

def postprocess(fname, field):
    """Postprocessing of the simulation."""
    if INPUTS['plotting']['solution']:
        fe.plot(field, title="Solution")
    if True in INPUTS['plotting'].values():
        fe.interactive()
    fe.File(fname+".pvd") << field

if __name__ == "__main__":
    ARGS = docopt(__doc__)
    if ARGS['-i']:
        INPUT_FILE = ARGS['-i']
    else:
        INPUT_FILE = 'simulation_inputs.json'
    with open(INPUT_FILE, 'r') as ifl:
        INPUTS = json.load(ifl)
    main()
