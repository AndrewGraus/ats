"""Richards water content evaluator: the standard form as a function of liquid saturation."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

""" Note that the density of rock could be either in mols or kg so this will
   only make sense if it's given in mols """

deps = [("porosity", "phi"),
        ("density_rock","nr"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("saturation_ice", "si"),
        ("molar_density_ice", "ni"),
        ("saturation_gas", "sg"),
        ("molar_density_gas", "ng")
        ]
params = []

import sympy
phi, nr, sl, nl, si, ni, sg, ng = sympy.var("phi,nr,sl,nl,si,ni,sg,ng")
expression =  nr*(1 - phi) + phi*(sl*nl + si*ni + sg*ng);

generate_evaluator("bulk_density", "ecosim_pk",
                   "bulk density", "bulk_density",
                   deps, params, expression=expression, doc=__doc__)
