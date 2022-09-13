import numpy as np
import sympy as sy
import sys

class Metric():
    """Create and manipulate a given covariant metric tensor.
    """

    def __init__(self, matrix, t, x, y=None, z=None):
        """Initialize the metric.

        parameter
        ---------
        matrix : numpy array (n x n)
            array containing the elements of the covariant metric tensor

        t, x, y, z : sympy symbols
            spacetime coordinates
        """

        shape = matrix.shape
        if len(shape) is not 2:
            sys.exit("Dimension error: the metric tensor is a 2-dim matrix")
        if shape[0] is not shape[1]:
            sys.exit("Dimension error: the metric tensor is a square matrix")

        self.dim = shape[0]                                                             # dimension of spacetime
        self.matrix = matrix                                                            # numpy-matrix of the metric tensor
        self.metric = sy.matrices.Matrix([col for col in matrix])                       # sympy-matrix of the metric tensor
        self.g = self.metric.det()                                                      # determinant of metric tensor
        self.inv_metric = self.metric.inv()                                             # inverse metric (sympy)
        self.t = t                                                                      # time coordinate
        self.x = x                                                                      # first space coordinate
        self.y = y                                                                      # second space coordinate
        self.z = z                                                                      # third space coordinate
        coords = [self.t, self.x, self.y, self.z]
        self.coords = [coord for coord in coords if coord is not None]                  # list of spacetime coordinates

        if len(self.coords) is not self.dim:
            sys.exit("dimension of metric and number of spacetime coordinates are not equal")

    def diff_metric(self, i, j, k):
        """Derivative of the component [i, j] of the metric tensor

        parameter
        ---------
        i, j : int
            indices of the specified element

        k : int
            derivative with respect to the spacetime coordinate with index k

        return
        ------
        deriv : sympy function
            derivative of the element
        """

        if i not in range(self.dim):
            sys.exit("index error in diff_metric: index i out of range")
        if j not in range(self.dim):
            sys.exit("index error in diff_metric: index j out of range")
        if k not in range(self.dim):
            sys.exit("index error in diff_metric: index k out of range")

        element = self.metric[i, j]
        deriv = sy.diff(element, self.coords[k])

        return deriv

    def christoffel_symbols(self, retC=True, simplify=True):
        """Calculation of all Christoffel symbols

        parameter
        ---------
        retC : bool
            if True the list of Christoffel symbols all_symbols will be returned

        return
        ------
        all_symbols : list (if retC=True)
            list containing all Christoffel symbols
        """

        # initialize a list with zero-matrices in the desired shape, such that all_symbols[k]_{ij} = Gamma^k_{ij}
        all_symbols = []
        for k in range(self.dim):
            symbol = sy.matrices.zeros(self.dim)
            all_symbols.append(symbol)

        # calculate the Christoffel symbols
        for k in range(self.dim):
            for i in range(self.dim):
                for j in range(self.dim):
                    for a in range(self.dim):
                        first_deriv = self.diff_metric(j, a, i)
                        second_deriv = self.diff_metric(i, a, j)
                        third_deriv = self.diff_metric(i, j, a)
                        all_symbols[k][i, j] += self.inv_metric[k, a] * (first_deriv + second_deriv - third_deriv) / 2

        # simplify if needed
        if simplify:
            for i in range(self.dim):
                all_symbols[i] = sy.simplify(all_symbols[i])

        # assign them to the class depending on self.dim
        if self.dim == 2:
            self.c0 = all_symbols[0]
            self.c1 = all_symbols[1]
            self.cs = [self.c0, self.c1]
        elif self.dim == 3:
            self.c0 = all_symbols[0]
            self.c1 = all_symbols[1]
            self.c2 = all_symbols[2]
            self.cs = [self.c0, self.c1, self.c2]
        else:
            self.c0 = all_symbols[0]
            self.c1 = all_symbols[1]
            self.c2 = all_symbols[2]
            self.c3 = all_symbols[3]
            self.cs = [self.c0, self.c1, self.c2, self.c3]

        if retC:
            return [symbol for symbol in all_symbols]

    def ricci_tensor(self, retR=True, simplify=True):
        """Calculation of the Ricci tensor

        parameter
        ---------
        retR : bool
            if True the function will return the Ricci tensor

        return
        ------
        riccitensor : matrix (if retR=True)
            Ricci tensor in sympy matrix format
        """

        # check whether Christoffel symbols exist or not
        try:
            self.c0
        except:
            self.christoffel_symbols(retC=False, simplify=False)

        # calculate Ricci tensor R_{ac}=R^d_{adc}
        riccitensor = sy.matrices.zeros(self.dim)
        for a in range(self.dim):
            for c in range(self.dim):
                for d in range(self.dim):
                    riccitensor[a, c] += sy.diff(self.cs[d][a, c], self.coords[d])
                    riccitensor[a, c] -= sy.diff(self.cs[d][a, d], self.coords[c])
                    for e in range(self.dim):
                        riccitensor[a, c] += self.cs[e][a, c] * self.cs[d][d, e]
                        riccitensor[a, c] -= self.cs[e][a, d] * self.cs[d][c, e]

        if simplify:
            riccitensor = sy.simplify(riccitensor)

        self.riccitensor = riccitensor

        if retR:
            return riccitensor

    def ricci_scalar(self, retR=True, simplify=True):
        """Calculation of the Ricci scalar

        parameter
        ---------
        retR : bool
            if True the function will return the Ricci scalar

        return
        ------
        ricciscalar : scalar (if retR=True)
            Ricci scalar in sympy matrix format
        """

        # Check whether Ricci tensor exists
        try:
            self.riccitensor
        except:
            self.ricci_tensor(retR=False, simplify=False)

        # calculate the Ricci scalar
        ricciscalar = 0
        contraction = self.inv_metric * self.riccitensor
        for i in range(self.dim):
            ricciscalar += contraction[i, i]

        if simplify:
            ricciscalar = sy.simplify(ricciscalar)

        self.ricciscalar = ricciscalar

        if retR:
            return ricciscalar

    def covariant_partial(self, retP=True):
        """A general co-/contravariant derivative has the form
        partial = a*partial_0 + b*partial_1 + c*partial_2 + d*partial_3
        for some constants a,b,c,d. For the partial derivative with covariant index we simply get
        partial_0 = a*partial_0 + b*partial_1 + c*partial_2 + d*partial_3 with a=1 and b=c=d=0. This simplification is not
        true for the partial derivative with contravariant index! Thats why the format of the partial is in this case maybe
        unnecessary difficult. The extra work will pay off in the next function 'contravariant_partial'.

        format
        ------
        [partial_mu] = [
            [1, 0, 0, 0],               <- mu=0
            [0, 1, 0, 0],               <- mu=1
            [0, 0, 1, 0],               <- mu=2
            [0, 0, 0, 1]                <- mu=3
        ]
             ^  ^  ^  ^
             |  |  |  |
             a  b  c  d

        parameter
        ---------
        retP : bool
            if True the function will return the covariant partial

        return
        ------
        covariantpartial : sympy matrix (if retP=True)
        """

        covariantpartial = sy.matrices.eye(self.dim)
        self.covariantpartial = covariantpartial

        if retP:
            return covariantpartial

    def contravariant_partial(self, retP=True):
        """Using the format from covariant_partial, the contravariant partial is simply given by the inverse metric itself!

        parameter
        ---------
        retP : bool
            if True the function will return the contravariant partial

        return
        ------
        contravariantpartial : sympy matrix (if retP=True)
        """

        contravariantpartial = self.inv_metric
        self.contravariantpartial = contravariantpartial

        if retP:
            return contravariantpartial

class MinkowskiMetric(Metric):
    def __init__(self):
        t, x, y, z = sy.symbols("t x y z", real=True)
        matrix = np.diag([1, -1, -1, -1])
        super(MinkowskiMetric, self).__init__(matrix=matrix, t=t, x=x, y=y, z=z)


class SchwarzschildMetric(Metric):
    def __init__(self):
        t, r, theta, phi, Rs = sy.symbols("t r theta phi Rs", real=True)
        matrix = np.array([
            [1-Rs/r, 0, 0, 0],
            [0, -1/(1-Rs/r), 0, 0],
            [0, 0, -r**2, 0],
            [0, 0, 0, -r**2*sy.sin(theta)**2]
        ])
        super(SchwarzschildMetric, self).__init__(matrix=matrix, t=t, x=r, y=theta, z=phi)
