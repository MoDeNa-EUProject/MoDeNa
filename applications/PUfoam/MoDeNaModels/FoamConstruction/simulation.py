#!/usr/bin/env python
"""Spatially three-dimensional model of heat/mass transfer.

Usage:
    simulation.py [-h | --help] [-i input_file] [--verbose]

Options:
    -h --help       Show this screen.
    -i input_file   Json file with inputs. Uses default file otherwise.
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
    system_matrix, field, bcs, cond = preprocess(fname)
    field = integrate(system_matrix, field, bcs)
    postprocess(fname, field, cond)
    print(
        term.yellow
        + "End."
        + term.normal
    )

def mesh_domain(domain):
    """Mesh computational domain using Gmsh."""
    call = sp.Popen(['gmsh', '-3', '-v', '3', domain])
    call.wait()

def convert_mesh(input_mesh, output_mesh):
    """Convert mesh to xml using dolfin-convert."""
    call = sp.Popen(['dolfin-convert', input_mesh, output_mesh])
    call.wait()

def bottomBC(x):
    """Bottom boundary condition."""
    return abs(x[2] - ZMAX) < fe.DOLFIN_EPS

def topBC(x):
    """Top boundary condition."""
    return abs(x[2] - ZMIN) < fe.DOLFIN_EPS

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


class Conductivity(fe.Expression):
    """
    Defines conductivity on computational domain. Maybe it should be rewritten
    as C++ code (see
    https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1005.html)
    """
    def __init__(self, subdomains, k_0, k_1, **kwargs):
        """Constructor."""
        self.subdomains = subdomains
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, cell):
        """
        Assigns value to the cell. Note that Gmsh indexes subdomains from 1.
        """
        if self.subdomains[cell.index] == 1:
            values[0] = self.k_0
        else:
            values[0] = self.k_1

def preprocess(fname):
    """Loads mesh, defines system of equations and prepares system matrix."""
    mesh = fe.Mesh(fname+".xml")
    if INPUTS['saving']['mesh']:
        fe.File(fname + "_mesh.pvd") << mesh
    # boundaries = fe.MeshFunction('size_t', mesh, fname+'_facet_region.xml')
    # if INPUTS['saving']['boundaries']:
    #     fe.File(fname+"_subdomains.pvd") << boundaries
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
        INPUTS['element_type'],
        INPUTS['element_degree'],
        constrained_domain=PeriodicDomain()
    )
    if ARGS['--verbose']:
        dofmap = fun_space.dofmap()
        print('Number of DOFs:', len(dofmap.dofs()))
    field = fe.TrialFunction(fun_space)
    test_func = fe.TestFunction(fun_space)
    cond = Conductivity(
        subdomains,
        fe.Constant(INPUTS['conductivity']['gas']),
        fe.Constant(INPUTS['conductivity']['solid']),
        degree=0
    )
    system_matrix = -cond*fe.inner(fe.grad(field), fe.grad(test_func))*fe.dx
    bctop = fe.Constant(INPUTS['boundary_conditions']['top'])
    bcbot = fe.Constant(INPUTS['boundary_conditions']['bottom'])
    # bcs = [
    #     fe.DirichletBC(fun_space, bctop, boundaries, 1),
    #     fe.DirichletBC(fun_space, bcbot, boundaries, 2)
    # ]
    bcs = [
        fe.DirichletBC(fun_space, bctop, topBC),
        fe.DirichletBC(fun_space, bcbot, bottomBC)
    ]
    field = fe.Function(fun_space)
    return system_matrix, field, bcs, cond

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

def postprocess(fname, field, cond):
    """Postprocessing of the simulation."""
    func_space = field.function_space()
    mesh = func_space.mesh()
    degree = func_space.ufl_element().degree()
    vec_func_space = fe.VectorFunctionSpace(
        mesh,
        INPUTS['element_type'],
        degree
    )
    flux = fe.project(-cond*fe.grad(field), vec_func_space)
    divergence = fe.project(-fe.div(cond*fe.grad(field)), func_space)
    flux_x, flux_y, flux_z = flux.split()
    av_flux = fe.assemble(flux_z*fe.dx)/((XMAX - XMIN)*(YMAX - YMIN))
    keff = av_flux*(ZMAX - ZMIN)/(
        INPUTS['boundary_conditions']['top']
        - INPUTS['boundary_conditions']['bottom']
    )
    print('Effective conductivity: {0}'.format(keff))
    with open(fname+"_keff.csv", 'w') as textfile:
        textfile.write('keff\n')
        textfile.write('{0}\n'.format(keff))
    fe.File(fname+"_solution.pvd") << field
    if INPUTS['saving']['flux']:
        fe.File(fname+"_flux.pvd") << flux
    if INPUTS['saving']['flux_divergence']:
        fe.File(fname+"_flux_divergence.pvd") << divergence
    if INPUTS['saving']['flux_components']:
        fe.File(fname+"_flux_x.pvd") << flux_x
        fe.File(fname+"_flux_y.pvd") << flux_y
        fe.File(fname+"_flux_z.pvd") << flux_z
    if INPUTS['plotting']['solution']:
        fe.plot(field, title="Solution")
    if INPUTS['plotting']['flux']:
        fe.plot(flux, title="Flux")
    if INPUTS['plotting']['flux_divergence']:
        fe.plot(divergence, title="Divergence")
    if INPUTS['plotting']['flux_components']:
        fe.plot(flux_x, title='x-component of flux (-kappa*grad(u))')
        fe.plot(flux_y, title='y-component of flux (-kappa*grad(u))')
        fe.plot(flux_z, title='z-component of flux (-kappa*grad(u))')
    if True in INPUTS['plotting'].values():
        fe.interactive()

if __name__ == "__main__":
    ARGS = docopt(__doc__)
    if ARGS['-i']:
        INPUT_FILE = ARGS['-i']
    else:
        INPUT_FILE = 'simulation_inputs.json'
    with open(INPUT_FILE, 'r') as ifl:
        INPUTS = json.load(ifl)
    main()
