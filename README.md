# phypy
A python library to automate certain calculations in physics
## Features
Note that this is only a list of the available features and not a documentation. For examples, see below.
Note furthermore, that calculations can be done symbolically with the help of sympy.
### General Relativity
The chosen sign convention is (+---).
Given the covariant metric tensor the following objects can be calculated automatically.
- contravariant metric tensor
- determinant of the covariant metric tensor
- Christoffel symbols
- Ricci tensor
- Ricci scalar
Predefined objects:
- Minkowski metric
- Schwarzschild metric
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
### Example: Christoffel symbols of the Schwarzschild metric
```python
metric = SchwarzschildMetric()
metric.christoffel_symbols(retC=False, simplify=True)
print(metric.c0)
```
By setting retC=False the function christoffel_symbols doent return the Christoffel symbols, but assign them to the class.
By simplify=True the sympy.simplify method will be used before assigning the Christoffel symbols to the class.
metric.c0 is the Christoffel symbol Gamma^0_{ij} in a matrix-format.
### Example: Define your own metric
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
