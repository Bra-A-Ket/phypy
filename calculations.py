from imports import *

def main():
    name = "phi"
    phi = RealScalarField4D(name=name, metric=FRWMetric())
    kg = phi.klein_gordon(retK=True, simplify=True, latex=False)
    print(kg)
if __name__ == "__main__":
    main()
