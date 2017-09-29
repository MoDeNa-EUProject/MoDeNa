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
from scipy import optimize
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
    # Load mesh and physical domains from file.
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
    # function space for temperature/concentration
    func_space = fe.FunctionSpace(
        mesh,
        INPUTS['element_type'],
        INPUTS['element_degree'],
        constrained_domain=PeriodicDomain()
    )
    # discontinuous function space for visualization
    dis_func_space = fe.FunctionSpace(
        mesh,
        'DG',
        INPUTS['element_degree'],
        constrained_domain=PeriodicDomain()
    )
    if ARGS['--verbose']:
        print('Number of cells:', mesh.num_cells())
        print('Number of faces:', mesh.num_faces())
        print('Number of edges:', mesh.num_edges())
        print('Number of vertices:', mesh.num_vertices())
        print('Number of DOFs:', len(func_space.dofmap().dofs()))
    # temperature/concentration field
    field = fe.TrialFunction(func_space)
    # test function
    test_func = fe.TestFunction(func_space)
    # function, which is equal to 1 everywhere
    unit_function = fe.Function(func_space)
    unit_function.assign(fe.Constant(1.0))
    # assign material properties to each domain
    if INPUTS['mode'] == 'conductivity':
        mat_prop = SubdomainConstant(
            subdomains,
            fe.Constant(INPUTS['conductivity']['gas']),
            fe.Constant(INPUTS['conductivity']['solid']),
            degree=0
        )
    elif INPUTS['mode'] == 'diffusivity':
        mat_prop = SubdomainConstant(
            subdomains,
            fe.Constant(INPUTS['diffusivity']['gas']),
            fe.Constant(INPUTS['diffusivity']['solid']) *
            fe.Constant(INPUTS['solubility'] *
                        gas_constant * INPUTS['temperature']),
            degree=0
        )
    # assign 1 to gas domain, and 0 to solid domain
    gas_content = SubdomainConstant(
        subdomains,
        fe.Constant(1.0),
        fe.Constant(0.0),
        degree=0
    )
    # define structure of foam over whole domain
    structure = fe.project(unit_function * gas_content, dis_func_space)
    # calculate porosity and wall thickness
    porosity = fe.assemble(structure *
                           fe.dx) / ((XMAX - XMIN) * (YMAX - YMIN) * (ZMAX - ZMIN))
    print('Porosity: {0}'.format(porosity))
    dwall = wall_thickness(
        porosity, INPUTS['morphology']['cell_size'],
        INPUTS['morphology']['strut_content'])
    print('Wall thickness: {0} m'.format(dwall))
    # calculate effective conductivity/diffusivity by analytical model
    if INPUTS['mode'] == 'conductivity':
        eff_prop = analytical_conductivity(
            INPUTS['conductivity']['gas'], INPUTS['conductivity']['solid'],
            porosity, INPUTS['morphology']['strut_content'])
        print('Analytical model: {0} W/(mK)'.format(eff_prop))
    elif INPUTS['mode'] == 'diffusivity':
        eff_prop = analytical_diffusivity(
            INPUTS['diffusivity']['solid'] *
            INPUTS['solubility'], INPUTS['solubility'],
            porosity, INPUTS['morphology']['cell_size'], dwall,
            INPUTS['temperature'], INPUTS['morphology']['enhancement_par'])
        print('Analytical model: {0} m^2/s'.format(eff_prop))
    # create system matrix
    system_matrix = -mat_prop * \
        fe.inner(fe.grad(field), fe.grad(test_func)) * fe.dx
    left_side, right_side = fe.lhs(system_matrix), fe.rhs(system_matrix)
    # define boundary conditions
    bcs = [
        fe.DirichletBC(func_space, fe.Constant(
            INPUTS['boundary_conditions']['top']), top_bc),
        fe.DirichletBC(func_space, fe.Constant(
            INPUTS['boundary_conditions']['bottom']), bottom_bc)
    ]
    # compute solution
    field = fe.Function(func_space)
    fe.solve(left_side == right_side, field, bcs)
    # output temperature/concentration at the boundaries
    if ARGS['--verbose']:
        print('Checking periodicity:')
        print('Value at XMIN:', field(XMIN, (YMIN + YMAX) / 3, (ZMIN + ZMAX) / 3))
        print('Value at XMAX:', field(XMAX, (YMIN + YMAX) / 3, (ZMIN + ZMAX) / 3))
        print('Value at YMIN:', field((XMIN + XMAX) / 3, YMIN, (ZMIN + ZMAX) / 3))
        print('Value at YMAX:', field((XMIN + XMAX) / 3, YMAX, (ZMIN + ZMAX) / 3))
        print('Value at ZMIN:', field((XMIN + XMAX) / 3, (YMIN + YMAX) / 3, ZMIN))
        print('Value at ZMAX:', field((XMIN + XMAX) / 3, (YMIN + YMAX) / 3, ZMAX))
    # calculate flux, and effective properties
    vec_func_space = fe.VectorFunctionSpace(
        mesh,
        INPUTS['element_type'],
        INPUTS['element_degree']
    )
    flux = fe.project(-mat_prop * fe.grad(field), vec_func_space)
    divergence = fe.project(-fe.div(mat_prop * fe.grad(field)), func_space)
    flux_x, flux_y, flux_z = flux.split()
    av_flux = fe.assemble(flux_z * fe.dx) / ((XMAX - XMIN) * (YMAX - YMIN))
    eff_prop = av_flux * (ZMAX - ZMIN) / (
        INPUTS['boundary_conditions']['top']
        - INPUTS['boundary_conditions']['bottom']
    )
    if INPUTS['mode'] == 'conductivity':
        print('Numerical model: {0} W/(mK)'.format(eff_prop))
    elif INPUTS['mode'] == 'diffusivity':
        print('Numerical model: {0} m^2/s'.format(eff_prop))
    # projection of concentration has to be in discontinuous function space
    if INPUTS['mode'] == 'diffusivity':
        sol_field = SubdomainConstant(
            subdomains,
            fe.Constant(1.0),
            fe.Constant(INPUTS['solubility'] *
                        gas_constant * INPUTS['temperature']),
            degree=0
        )
        field = fe.project(field * sol_field, dis_func_space)
    # save results
    with open(fname + "_eff_prop.csv", 'w') as textfile:
        textfile.write('eff_prop\n')
        textfile.write('{0}\n'.format(eff_prop))
    fe.File(fname + "_solution.pvd") << field
    fe.File(fname + "_structure.pvd") << structure
    if INPUTS['saving']['flux']:
        fe.File(fname + "_flux.pvd") << flux
    if INPUTS['saving']['flux_divergence']:
        fe.File(fname + "_flux_divergence.pvd") << divergence
    if INPUTS['saving']['flux_components']:
        fe.File(fname + "_flux_x.pvd") << flux_x
        fe.File(fname + "_flux_y.pvd") << flux_y
        fe.File(fname + "_flux_z.pvd") << flux_z
    # plot results
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


def analytical_diffusivity(perm, sol, por, dcell, dwall, temp, enh):
    """Effective diffusivity."""
    return enh * dcell / dwall * perm * (por / (gas_constant * temp)
                                         + sol * (1 - por))**(-1)


def wall_thickness(por, dcell, fstrut):
    """Determines wall thickness. Placido."""
    sol = optimize.root(wall_thickness_root, [
        1e-6, 1e-6], args=(por, dcell, fstrut))
    return sol.x[0]


def wall_thickness_root(x, *args):
    """Root function for the determination of wall thickness. Placido."""
    dwall = x[0]
    dstrut = x[1]
    por = args[0]
    dcell = args[1]
    fstrut = args[2]
    vcell = 0.348 * dcell**3
    vstruts = 2.8 * dstrut**2 * dcell - 3.93 * dstrut**3
    vwalls = (1.3143 * dcell**2 - 7.367 * dstrut *
              dcell + 10.323 * dstrut**2) * dwall
    return [vstruts / (vstruts + vwalls) - fstrut,
            ((1 - por) * 0.348 * dcell**3 - 2.8 * dstrut**2 * dcell
             + 3.93 * dstrut**3) / (1.3143 * dcell**2 - 7.367 * dstrut * dcell
                                    + 10.323 * dstrut**2) - dwall]


def wall_thickness2(por, dcell, fstrut):
    """Determines wall thickness. Kaemmerlen."""
    sol = optimize.root(wall_thickness_root2, [
        1e-6, 1e-6], args=(por, dcell, fstrut))
    return sol.x[0]


def wall_thickness_root2(x, *args):
    """Root function for the determination of wall thickness. Kaemmerlen."""
    dwall = x[0]
    dstrut = x[1]
    por = args[0]
    dcell = args[1]
    fstrut = args[2]
    vcell = 0.349 * dcell**3
    vstruts = 2.805 * dstrut**2 * dcell
    vwalls = (1.317 * dcell**2 - 13.4284 * dstrut *
              dcell + 34.2375 * dstrut**2) * dwall + (4.639 * dcell
                                                      - 17.976 * dstrut) * dwall**2
    return [vstruts / (vstruts + vwalls) - fstrut,
            vwalls + vstruts - (1 - por) * vcell]


if __name__ == "__main__":
    ARGS = docopt(__doc__)
    if ARGS['-i']:
        INPUT_FILE = ARGS['-i']
    else:
        INPUT_FILE = 'simulation_inputs.json'
    with open(INPUT_FILE, 'r') as ifl:
        INPUTS = json.load(ifl)
    main()
