from imports import *

def main():
    names = ["phi_1", "phi_2", "phi_3", "phi_3"]
    fields = [RealScalarField4D(name=name) for name in names]
    WickContraction(fields=fields)
if __name__ == "__main__":
    main()
