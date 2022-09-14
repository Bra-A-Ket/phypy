from imports import *

def main():
    name = "phi"
    t, r, theta, phi = sy.symbols("t r theta phi", real=True)
    phi = RealScalarField4D(name=name, metric=FRWMetric())
    kg = phi.klein_gordon(retK=True, simplify=True, latex=False)
    print(kg)
if __name__ == "__main__":
    main()
