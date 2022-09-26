# phypy
A python library to automate certain calculations in physics
## Features
Note that this is only a list of the available features and not a documentation. Some examples can be found below. For
additional information visit the [Wiki](https://github.com/Bra-A-Ket/phypy/wiki).

Note furthermore, that calculations can be done symbolically with the help of sympy.
### General Relativity
The chosen sign convention is (+---).
Given the covariant metric tensor the following objects can be calculated automatically.
- contravariant metric tensor
- determinant of the covariant metric tensor
- Christoffel symbols
- Ricci tensor
- Ricci scalar
- co- and contravariant partial derivative (partial_mu / partial^mu)

Predefined objects:
- Minkowski metric
- Schwarzschild metric
- Friedmann-Robertson-Walker metric
### Quantum Field Theory
If no metric is specified the Minkowski metric is assumed.
- Minimal coupled Klein-Gordon equation in the background of a given GR metric
- Wick contractions
## Requirements
In order to make phypy work, one has to install some external packages.
- numpy for obvious reasons
```console
python -m pip install numpy
```
- sympy for symbolic calculations
```console
python -m pip install sympy
```
## Usage
Feel free to use any python file, e.g. calculations.py, and make sure that you have imported everything:
```python
from imports import *
```
### Example 1: Christoffel symbols of the Schwarzschild metric
```python
metric = SchwarzschildMetric()
metric.christoffel_symbols(retC=False, simplify=True)
print(metric.c0)
```
By setting retC=False the function christoffel_symbols doent return the Christoffel symbols, but assign them to the class.
By simplify=True the sympy.simplify method will be used before assigning the Christoffel symbols to the class.
metric.c0 is the Christoffel symbol Gamma^0_{ij} in a matrix-format.
### Example 2: Define your own metric
```python
t, x, y, z = sy.symbols("t x y z", real=True)
matrix = np.array([
    [t, 0, 0, 0],
    [0, -x, 0, 0],
    [0, 0, -y, 0],
    [0, 0, 0, -z**2]
])
metric = Metric(matrix=matrix, t=t, x=x, y=y, z=z)
metric.ricci_scalar(retR=False, simplify=True)
print(metric.ricciscalar)
```
### Example 3: Ricci scalar of FRW metric
```python
metric = FRWMetric()
ricciscalar = metric.ricci_scalar(retR=True, simplify=True)
print(ricciscalar)
```
### Example 4: Klein-Gordon equation in FRW metric
```python
name = "phi"
phi = RealScalarField4D(name=name, metric=FRWMetric())
kg = phi.klein_gordon(retK=True, simplify=True, latex=False)
print(kg)
```
### Example 5: Wick contractions
Find all Wick contractions of <0|T phi_1 phi_2 phi_3 phi_3|0>.
```python
names = ["phi_1", "phi_2", "phi_3", "phi_3"]
fields = [RealScalarField4D(name=name) for name in names]
WickContraction(fields=fields)
```
The output looks like this:
```console
<0|T[1, 2, 3, 3]|0> =

1 x [[1, 2], [3, 3]] +
2 x [[1, 3], [2, 3]]
```
This should be read as: <0|T phi_1 phi_2 phi_3 phi_3|0> = 1x<0|T phi_1 phi_2|0><0|T phi_3 phi_3|0> + 2x<0|T phi_1 phi_3>
<0|T phi_2 phi_3|0>.
