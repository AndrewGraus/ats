"""Richards water content evaluator: the standard form as a function of liquid saturation."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("porosity", "phi"),
        ("saturation_liquid", "sl"),
        ("molar_density_liquid", "nl"),
        ("cell_volume", "cv")
        ]

params = [("alpha", "double", "van genuchten alpha"),
("m", "double", "van genuchten m"),
("n", "double", "van genuchten n"),
("sr", "double", "residual saturation")
]

import sympy
phi, sl, nl, cv = sympy.var("phi,sl,nl,cv")
alpha_, m_, n_, sr_ = sympy.var("alpha_,m_,n_,sr_")
expression = -1/(alpha_**n_) * (1 - ((sl - sr_)/(cv*nl*phi - sr_))**(1/m_))**(-n_)

generate_evaluator("matric_pressure", "flow",
                   "matric pressure", "matric_pressure",
                   deps, params, expression=expression, doc=__doc__)
