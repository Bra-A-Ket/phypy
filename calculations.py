from imports import *

def main():

    t, x, y, z = sy.symbols("t x y z", real=True)
    N = sy.Function("N")(t, x, y, z)
    G11 = sy.Function("G_{11}")(t, x, y, z)
    G12 = sy.Function("G_{12}")(t, x, y, z)
    G13 = sy.Function("G_{13}")(t, x, y, z)
    G22 = sy.Function("G_{22}")(t, x, y, z)
    G23 = sy.Function("G_{23}")(t, x, y, z)
    G33 = sy.Function("G_{33}")(t, x, y, z)
    matrix = np.array([
        [N**2, 0 , 0 , 0],
        [0, G11, G12, G13],
        [0, G12, G22, G23],
        [0, G13, G23, G33]
    ])
    metric = Metric(matrix=matrix, t=t, x=x, y=y, z=z)

    name = "phi"
    phi = RealScalarField4D(name=name)
    phi.klein_gordon(retK=False, simplify=True, latex=False, metric=metric)
if __name__ == "__main__":
    main()
