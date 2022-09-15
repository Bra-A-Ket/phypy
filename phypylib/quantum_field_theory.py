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
    """Real scalar field in 4 spacetime dimensions.

    parameter
    ---------
    name : string
        name of the field

    m : sympy symbol
        symbol for the mass term. Can also be set m=0

    metric : spacetime metric from general_relativity.Metric()
        specified metric of spacetime

    # TODO: 3- and 2-dim realscalarfield using this class -> edit MinkowskiMetric()
    """

    def __init__(self, name, m=sy.symbols("m", real=True), metric=MinkowskiMetric()):
        try:
            isinstance(name, str)
        except:
            sys.exit("RealScalarField4D: name variable is a string")

        self.m = m                                                                      # mass term
        self.fieldtype = "4-dim real scalar field"                                      # field type
        self.metric = metric                                                            # metric
        t = self.metric.t
        x = self.metric.x
        y = self.metric.y
        z = self.metric.z
        field = sy.Function(name, real=True)(t, x, y, z)
        super(RealScalarField4D, self).__init__(field=field, x=x, y=y, z=z, t=t)

    def gr_dalembert_operator(self, retG=True, simplify=True, latex=False):
        """Calculates the d'Alembert operator acting on the field Nabla_mu*Nabla^mu*field for a given metric.

        parameter
        ---------
        retG : bool
            if True the result is returned

        simplify : bool
            if True the sympy.simplify function is used

        latex : bool
            if True the result will be printed in latex format

        return
        ------
        dalembert
        """

        covariantpartial = self.metric.covariant_partial(retC=True)
        contravariantpartial = self.metric.contravariant_partial(retC=True)

        # dalembert acting on phi
        dalembert = 0
        for i in range(self.dim):
            for j in range(self.dim):
                dalembert += sy.diff(contravariantpartial[i, j]*sy.diff(self.field, self.coords[j]), self.coords[i])
        if simplify:
            dalembert = sy.simplify(dalembert)

        # correction term from Christoffel symbols
        self.metric.christoffel_symbols(retC=False)
        for i in range(self.dim):
            for j in range(self.dim):
                dalembert += self.metric.cs[i][i, j] * contravariantpartial[i, j] * sy.diff(self.field, self.coords[j])

        if simplify:
            dalembert = sy.simplify(dalembert)
        self.dalembert = dalembert

        if latex:
            sy.print_latex(dalembert)

        if retG:
            return dalembert

    def klein_gordon(self, retK=True, simplify=True, latex=False):
        """Calculate the Klein-Gordon equation for a massive real scalar field in 4D using the given metric by initializing
        the field.

        parameter
        ---------
        retK : bool
            if True the result is returned

        simplify : bool
            if True the sympy.simplify function is used

        latex : bool
            if True the result will be printed in latex format

        return
        ------
        kleingordon
            how to read: 0 = kleingordon
        """
        self.gr_dalembert_operator(retG=False, simplify=simplify, latex=False)

        kleingordon = self.dalembert + self.m**2*self.field
        self.kleingordon = kleingordon

        if latex:
            sy.print_latex(kleingordon)

        if retK:
            return kleingordon
