#!/usr/bin/env python
"""Determination of initial cell gas composition.

Usage:
    cell_gas.py [-h | --help] [-i input_file] [--verbose]

Options:
    -h --help       Show this screen.
    -i input_file   Json file with inputs. Uses default file otherwise.
    --verbose       Print more information.

@author: Pavel Ferkl
"""
from __future__ import division, print_function
import csv
import json
import os
from docopt import docopt
from scipy.constants import gas_constant
from numpy import exp, linspace
from modena import SurrogateModel


def saturated_pressure(name, temp):
    """
    Saturated pressure of gas. Implemented for cyclopentane, carbon dioxide,
    Opteon 1100, and Solstice.
    http://webbook.nist.gov/cgi/cbook.cgi?ID=C287923&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
    http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
    Product sheets are used as data source for Opteon and Solstice.
    """
    if name == 'CyP':
        if temp < 288:
            par = [
                4.24714,
                1235.305,
                -30.666
            ]
        else:
            par = [
                4.00288,
                1119.208,
                -42.412
            ]
    elif name == 'CO2':
        """Validity 154-196 K, but function should ensure that CO2 does not
        condense."""
        par = [
            6.81228,
            1301.679,
            -3.494
        ]
    elif name == 'Opt':
        par = [
            4.72048976,
            1493.15451888,
            17.71674892
        ]
    elif name == 'Sol':
        par = [
            4.18988677,
            1071.53370906,
            -35.18925819
        ]
    else:
        raise Exception('Saturated pressure for {} undefined'.format(name))
    return 1e5 * 10**(par[0] - (par[1] / (temp + par[2])))


def surface_tension(name, temp):
    """
    Surface tension of liquid. Implemented for cyclopentane and carbon dioxide.
    Uses cyclopentane data otherwise.
    http://www.ddbst.com/en/EED/PCP/SFT_C1050.php
    http://linkinghub.elsevier.com/retrieve/pii/S0021961405002314
    """
    if name == 'CO2':
        # Data for 217-302 K. Extrapolated for higher temperatures.
        if temp < 302.136780079:
            par = [
                1.05646453e+02,
                - 5.60265519e-01,
                6.97039246e-04
            ]
        else:
            par = [
                0, 0, 0
            ]
    else:
        # Data for 283-313 K. Extrapolated for higher temperatures.
        if temp < 408.482827604:
            par = [
                2.21308934e+01,
                1.45524444e-01,
                - 4.88888889e-04
            ]
        else:
            par = [
                0, 0, 0
            ]
    return 1e-3 * (par[0] + par[1] * temp + par[2] * temp**2)


def kelvin_effect(pres, surft, temp, mw_ba, dcell):
    """Saturated vapour pressure changes with curvature. It is assumed that
    cells are spherical. https://en.wikipedia.org/wiki/Kelvin_equation. This
    effect can be generally neglected when curvature radius is smaller than 100
    nm."""
    volm = mw_ba/1e3 # approximation: using density 1000 kg/m3
    return pres*exp(-4*surft*volm/(dcell*gas_constant*temp))


def solubility_in_polymer(name, temp):
    """solubility in polymer"""
    model = SurrogateModel.load('Solubility[A=' + name + ',B=2]')
    inputs = {'T': temp, 'xl1': 1e-3, 'xl2': 0.999}
    outputs = model.callModel(inputs)
    henry = outputs['H']
    return henry


def initial_pressure(name, w_ba_ini, temp, args):
    """Calculates gaseous, dissolved, and condensed weight fractions."""
    por, mw_ba, m_foam, m_pol, volume, dcell = args
    pres_ba = w_ba_ini * m_foam / (  # partial pressure of blowing agent
        mw_ba[name] * volume * por / gas_constant / temp
        + solubility_in_polymer(name, temp) * m_pol / 1e5
    )
    pres_ba_sat = saturated_pressure(name, temp)  # vapor pressure
    surft = surface_tension(name, temp) # surface tension
    pres_ba_sat = kelvin_effect(pres_ba_sat, surft, temp, mw_ba[name], dcell)
    # pressure of blowing agent cannot be higher than vapor pressure
    if pres_ba > pres_ba_sat:
        pres_ba = pres_ba_sat
    if ARGS['--verbose']:
        print('Temperature {0:.1f} K'.format(temp))
        print('Partial pressure of blowing agent: {0:.0f} Pa'.format(pres_ba))
        print('Vapor pressure: {0:.0f} Pa'.format(pres_ba_sat))
    w_ba_g = mw_ba[name] * pres_ba * volume * \
        por / gas_constant / temp / m_foam
    w_ba_d = solubility_in_polymer(
        name, temp) * m_pol * pres_ba / 1e5 / m_foam
    w_ba_c = w_ba_ini - w_ba_g - w_ba_d
    if w_ba_c < 1e-8:  # nicer output
        w_ba_c = 0.0
    if ARGS['--verbose']:
        print('Weight fraction of BA in gas phase: {0:.3g} g/g'.format(w_ba_g))
        print('Weight fraction of BA dissolved in polymer: {0:.3g} g/g'.format(
            w_ba_d
        ))
        print('Weight fraction of BA condensed: {0:.3g} g/g'.format(w_ba_c))
    return pres_ba, w_ba_g, w_ba_d, w_ba_c


def main():
    """Main function. Executed from command line."""
    print('Calculating gas compositions.')
    resf = 'results/cell_gas' # results folder
    if not os.path.isdir(resf):
        os.makedirs(resf)
    sizex = 0.03  # sample size
    sizey = 0.02  # sample size
    sizez = 0.02  # sample size
    volume = sizex**3  # sample volume
    volume = sizex * sizey * sizez  # sample volume
    # polymer density
    rhop = INPUTS['polymer_density']
    # molecular weight
    mw_ba = INPUTS['molar_mass']
    # foam density
    rhof = INPUTS['foam_density']
    # cell size for Kelvin effect on saturated vapour pressure
    dcell = INPUTS['cell_size']
    # initial weight fraction of BA
    w_ba_ini = INPUTS['initial_weight_fraction']
    names = w_ba_ini.keys()
    if 'H2O' in w_ba_ini:
        if 'CO2' in w_ba_ini:
            print("WARNING: H2O and CO2 are both in initial_weight_fraction.",
                  "We will sum these contributions.")
        else:
            w_ba_ini['CO2'] = 0
        w_ba_ini['CO2'] += w_ba_ini['H2O'] * mw_ba['CO2'] / mw_ba['H2O']
        names.append('CO2')
        names.remove('H2O')
    temps = linspace(
        INPUTS['temperature']['min'],
        INPUTS['temperature']['max'],
        INPUTS['temperature']['points']
    )
    por = 1 - rhof / rhop  # porosity
    m_foam = rhof * volume  # foam sample weight
    m_pol = m_foam * (1 - sum(w_ba_ini.values()))  # weight of polymer
    if ARGS['--verbose']:
        print('Foam weight {0:.3f} g'.format(m_foam * 1e3))
    args = [por, mw_ba, m_foam, m_pol, volume, dcell]
    for name in names:
        with open(os.path.join(resf, 'cell_gas_{0}.csv'.format(name)),
                  'w') as csvfile:
            fieldnames = ['temp', 'pres_ba', 'w_ba_g', 'w_ba_d', 'w_ba_c']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for temp in temps:
                pres_ba, w_ba_g, w_ba_d, w_ba_c = initial_pressure(
                    name, w_ba_ini[name], temp, args)
                writer.writerow(
                    {'temp': temp, 'pres_ba': pres_ba, 'w_ba_g': w_ba_g,
                     'w_ba_d': w_ba_d, 'w_ba_c': w_ba_c})
    print('End.')


if __name__ == "__main__":
    ARGS = docopt(__doc__)
    if ARGS['-i']:
        INPUT_FILE = ARGS['-i']
    else:
        INPUT_FILE = 'inputs/cell_gas.json'
    with open(INPUT_FILE, 'r') as ifl:
        INPUTS = json.load(ifl)
    main()
