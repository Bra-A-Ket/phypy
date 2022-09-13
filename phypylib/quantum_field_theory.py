import numpy as np
import sympy as sy
import sys
from phypylib.general_relativity import *

class Field():
    def __init__(self, field, x, y=None, z=None, t=None):
        self.field = field                                                              # field function
        self.t = t                                                                      # time coordinate
        self.x = x                                                                      # first space coordinate
        self.y = y                                                                      # second space coordinate
        self.z = z                                                                      # third space coordinate
        coords = [self.t, self.x, self.y, self.z]
        self.coords = [coord for coord in coords if coord is not None]                  # list of spacetime coordinates
        self.dim = len(self.coords)


class RealScalarField4D(Field):
    def __init__(self, name, m=sy.symbols("m", real=True)):
        try:
            isinstance(name, str)
        except:
            sys.exit("RealScalarField4D: name variable is a string")

        self.m = m
        self.fieldtype = "4-dim real scalar field"
        t, x, y, z = sy.symbols("t x y z", real=True)
        field = sy.Function(name, real=True)(t, x, y, z)
        super(RealScalarField4D, self).__init__(field=field, x=x, y=y, z=z, t=t)

    def klein_gordon(self, metric=MinkowskiMetric(), retK=True, simplify=True, latex=False):
        covariantpartial = metric.covariant_partial(retC=True)
        contravariantpartial = metric.contravariant_partial(retC=True)

        # dynamic part {nabla_mu*nabla^mu field}
        dalembert = 0
        for i in range(self.dim):
            for j in range(self.dim):
                dalembert += sy.diff(contravariantpartial[i, j]*sy.diff(self.field, self.coords[j]), self.coords[i])
        if simplify:
            dalembert = sy.simplify(dalembert)
        self.dalembert = dalembert

        kleingordon = self.dalembert + self.m**2*self.field
        metric.christoffel_symbols(retC=False)
        for i in range(self.dim):
            for j in range(self.dim):
                kleingordon += metric.cs[i][i, j] * contravariantpartial[i, j] * sy.diff(self.field, self.coords[j])

        if simplify:
            kleingordon = sy.simplify(kleingordon)
        self.kleingordon = kleingordon

        if latex:
            sy.print_latex(kleingordon)

        if retK:
            return kleingordon
