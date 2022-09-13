from imports import *

def main():
    """
    theta, phi = sy.symbols("theta, phi", real=True)
    R = sy.symbols("R", real=True)
    matrix = np.array([
        [R**2, 0],
        [0, R**2*sy.sin(theta)**2]
    ])
    metric = Metric(matrix=matrix, t=theta, x=phi)
    metric.ricci_scalar(retR=False)
    print(metric.ricciscalar)
    """


    schwarzschild_metric = SchwarzschildMetric()
    schwarzschild_metric.christoffel_symbols(retC=False, simplify=False)
    print(schwarzschild_metric.c0)
    """
    t, x, y, z = sy.symbols("t x y z", real=True)
    N = sy.Function("N")(t, x, y, z)
    G11 = sy.Function("G11")(t, x, y, z)
    G12 = sy.Function("G12")(t, x, y, z)
    G13 = sy.Function("G13")(t, x, y, z)
    G22 = sy.Function("G22")(t, x, y, z)
    G23 = sy.Function("G23")(t, x, y, z)
    G33 = sy.Function("G33")(t, x, y, z)
    matrix = np.array([
        [N**2, 0 , 0 , 0],
        [0, G11, G12, G13],
        [0, G12, G22, G23],
        [0, G13, G23, G33]
    ])
    metric = Metric(matrix=matrix, t=t, x=x, y=y, z=z)
    metric.christoffel_symbols(retC=False)
    print(metric.c0)
    """

if __name__ == "__main__":
    main()
