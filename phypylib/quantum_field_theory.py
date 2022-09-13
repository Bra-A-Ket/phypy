import numpy as np
import sympy as sy
import sys

class Field():
    def __init__(self, name, x, y=None, z=None, t=None):
        self.name = name                                                                # field variable name
        self.t = t                                                                      # time coordinate
        self.x = x                                                                      # first space coordinate
        self.y = y                                                                      # second space coordinate
        self.z = z                                                                      # third space coordinate
        coords = [self.t, self.x, self.y, self.z]
        self.coords = [coord for coord in coords if coord is not None]                  # list of spacetime coordinates


class RealScalarField4D(Field):
    def __init__(self, name):
        name = sy.symbols(name, real=True)
        t, x, y, z = sy.symbols("t x y z", real=True)
        super(RealScalarField4D, self).__init__(name=name, x=x, y=y, z=z, t=t)
