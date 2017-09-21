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
import fenics as fe
from blessings import Terminal
from docopt import docopt
from scipy.constants import gas_constant
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
    system_matrix, field, bcs, diff, subdomains = preprocess(fname)
    field = integrate(system_matrix, field, bcs)
    postprocess(fname, field, diff, subdomains)
    print(
        term.yellow
        + "End."
        + term.normal
    )


def bottom_bc(x):
    """Bottom boundary condition."""
    return abs(x[2] - ZMAX) < fe.DOLFIN_EPS


def top_bc(x):
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


class SubdomainConstant(fe.Expression):
    """
    Defines constant on computational domain respecting the subdomains. Maybe it
    should be rewritten as C++ code (see
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


def analytical_conductivity(k_gas, k_sol, por, f_strut):
    """Effective conductivity according to Ahern et al. (2005)."""
    fun = 2 * (1 - f_strut) / 3 * (1 + k_gas / (2 * k_sol)) + \
        f_strut / 3 * (1 + 4 * k_gas / (k_gas + k_sol))
    return (k_gas * por + k_sol * fun * (1 - por)) / (por + (1 - por) * fun)


def analytical_diffusivity(perm, sol, por, dcell, dwall, temp, geom):
    """Effective diffusivity."""
    rgas = 8.314
    return geom * dcell / dwall * perm * (por / (rgas * temp)
                                          + sol * (1 - por))**(-1)
    # return geom * dcell / dwall * perm * (por / (rgas * temp))**(-1)


def preprocess(fname):
    """Loads mesh, defines system of equations and prepares system matrix."""
    mesh = fe.Mesh(fname + ".xml")
    if INPUTS['saving']['mesh']:
        fe.File(fname + "_mesh.pvd") << mesh
    if INPUTS['plotting']['mesh']:
        fe.plot(mesh, title='Mesh')
    subdomains = fe.MeshFunction(
        'size_t', mesh, fname + '_physical_region.xml')
    if INPUTS['saving']['subdomains']:
        fe.File(fname + "_subdomains.pvd") << subdomains
    if INPUTS['plotting']['subdomains']:
        fe.plot(subdomains, title='Subdomains')
    # Boundaries cannot be exported from gmsh for some foams.
    # Possible bug in gmsh? Boundaries are defined through functions here instead.
    # boundaries = fe.MeshFunction('size_t', mesh, fname+'_facet_region.xml')
    # if INPUTS['saving']['boundaries']:
    #     fe.File(fname+"_subdomains.pvd") << boundaries
    # if INPUTS['plotting']['boundaries']:
    #     fe.plot(boundaries, title='Boundaries')
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
    if INPUTS['mode'] == 'conductivity':
        diff = SubdomainConstant(
            subdomains,
            fe.Constant(INPUTS['conductivity']['gas']),
            fe.Constant(INPUTS['conductivity']['solid']),
            degree=0
        )
    elif INPUTS['mode'] == 'diffusivity':
        diff = SubdomainConstant(
            subdomains,
            fe.Constant(INPUTS['diffusivity']['gas']),
            fe.Constant(INPUTS['diffusivity']['solid']) *
            fe.Constant(INPUTS['solubility'] *
                        gas_constant * INPUTS['temperature']),
            degree=0
        )
    unit_function = fe.Function(fun_space)
    unit_function.assign(fe.Constant(1.0))
    gas_content = SubdomainConstant(
        subdomains,
        fe.Constant(1.0),
        fe.Constant(0.0),
        degree=0
    )
    porosity = fe.assemble(gas_content * unit_function *
                           fe.dx) / (XMAX * YMAX * ZMAX)
    print('Porosity: {0}'.format(porosity))
    if INPUTS['analytical_model']:
        if INPUTS['mode'] == 'conductivity':
            keff = analytical_conductivity(
                INPUTS['conductivity']['gas'], INPUTS['conductivity']['solid'],
                porosity, 0.0)
            print('Analytical model: {0} W/(mK)'.format(keff))
        elif INPUTS['mode'] == 'diffusivity':
            deff = analytical_diffusivity(
                INPUTS['diffusivity']['solid'] *
                INPUTS['solubility'], INPUTS['solubility'],
                porosity, INPUTS['morphology']['cell_size'],
                INPUTS['morphology']['wall_thickness'],
                INPUTS['temperature'], 1.0)
            print('Analytical model: {0} m^2/s'.format(deff))
    system_matrix = -diff * \
        fe.inner(fe.grad(field), fe.grad(test_func)) * fe.dx
    bctop = fe.Constant(INPUTS['boundary_conditions']['top'])
    bcbot = fe.Constant(INPUTS['boundary_conditions']['bottom'])
    # bcs = [
    #     fe.DirichletBC(fun_space, bctop, boundaries, 1),
    #     fe.DirichletBC(fun_space, bcbot, boundaries, 2)
    # ]
    bcs = [
        fe.DirichletBC(fun_space, bctop, top_bc),
        fe.DirichletBC(fun_space, bcbot, bottom_bc)
    ]
    field = fe.Function(fun_space)
    return system_matrix, field, bcs, diff, subdomains


def integrate(system_matrix, field, bcs):
    """Integrates the equation"""
    left_side, right_side = fe.lhs(system_matrix), fe.rhs(system_matrix)
    fe.solve(left_side == right_side, field, bcs)
    if ARGS['--verbose']:
        print('Checking periodicity:')
        print('Value at XMIN:', field(XMIN, (YMIN + YMAX) / 3, (ZMIN + ZMAX) / 3))
        print('Value at XMAX:', field(XMAX, (YMIN + YMAX) / 3, (ZMIN + ZMAX) / 3))
        print('Value at YMIN:', field((XMIN + XMAX) / 3, YMIN, (ZMIN + ZMAX) / 3))
        print('Value at YMAX:', field((XMIN + XMAX) / 3, YMAX, (ZMIN + ZMAX) / 3))
        print('Value at ZMIN:', field((XMIN + XMAX) / 3, (YMIN + YMAX) / 3, ZMIN))
        print('Value at ZMAX:', field((XMIN + XMAX) / 3, (YMIN + YMAX) / 3, ZMAX))
    return field


def postprocess(fname, field, diff, subdomains):
    """Postprocessing of the simulation."""
    func_space = field.function_space()
    mesh = func_space.mesh()
    degree = func_space.ufl_element().degree()
    vec_func_space = fe.VectorFunctionSpace(
        mesh,
        INPUTS['element_type'],
        degree
    )
    flux = fe.project(-diff * fe.grad(field), vec_func_space)
    divergence = fe.project(-fe.div(diff * fe.grad(field)), func_space)
    flux_x, flux_y, flux_z = flux.split()
    av_flux = fe.assemble(flux_z * fe.dx) / ((XMAX - XMIN) * (YMAX - YMIN))
    keff = av_flux * (ZMAX - ZMIN) / (
        INPUTS['boundary_conditions']['top']
        - INPUTS['boundary_conditions']['bottom']
    )
    if INPUTS['mode'] == 'conductivity':
        print('Numerical model: {0} W/(mK)'.format(keff))
    elif INPUTS['mode'] == 'diffusivity':
        print('Numerical model: {0} m^2/s'.format(keff))
    # projection has to be to discontinuous function space
    if INPUTS['mode'] == 'diffusivity':
        sol_field = SubdomainConstant(
            subdomains,
            fe.Constant(1.0),
            fe.Constant(INPUTS['solubility'] *
                        gas_constant * INPUTS['temperature']),
            degree=0
        )
        dis_func_space = fe.FunctionSpace(
            mesh,
            'DG',
            INPUTS['element_degree'],
            constrained_domain=PeriodicDomain()
        )
        field = fe.project(field * sol_field, dis_func_space)
    with open(fname + "_keff.csv", 'w') as textfile:
        textfile.write('keff\n')
        textfile.write('{0}\n'.format(keff))
    fe.File(fname + "_solution.pvd") << field
    if INPUTS['saving']['flux']:
        fe.File(fname + "_flux.pvd") << flux
    if INPUTS['saving']['flux_divergence']:
        fe.File(fname + "_flux_divergence.pvd") << divergence
    if INPUTS['saving']['flux_components']:
        fe.File(fname + "_flux_x.pvd") << flux_x
        fe.File(fname + "_flux_y.pvd") << flux_y
        fe.File(fname + "_flux_z.pvd") << flux_z
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
