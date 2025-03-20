# Coordinate differentiation for Mathematica
This package/library provides tools for performing differentiation and index contractions on tensor structures in $\mathbb{R}^d$. The structures are assumed to respect $O(d)$ symmetry and therefore depend only on coordinate vectors of different points and the distances between them. The Python library is built on top of SymPy.
## Installation

### Python
Clone this repository using
```console
git clone https://github.com/AlexanderGTumanov/RdDiff.git
```
Navigate to the ``RdDiff`` directory:
```console
cd RdDiff
```
install the package:
```console
pip install -e .
```
Now you can import it in Python:
```console
import RdDiff as rdd
from RdDiff import x
from RdDiff import xx
from RdDiff import d
from RdDiff import delta
```
### Mathematica
Inside the repository, you'll find the ``mathematica`` folder with ``RdDiff.m``.
In Mathematica:
```mathematica
Get["/path/to/RdDiff/mathematica/RdDiff.m"]
```
Replace ``/path/to`` with the actual path where you cloned the repository.
## Usage
This package is designed to perform various operations on d-dimensional functions, vectors, and tensors that depend on multiple coordinate vectors. The metric signature is irrelevant, as the code does not rely on it. The package's notation can be summarized in the following table,

| Object                                | Python notations    |  Mathematica notations  |
| --------                              | -------             | -------                 |
| $x^i_a$                               | ``x(a, i)``         | ``x[a,i]``              |
| $x_a^2 = x_a^i\ x_{ai}$               | ``xx(a)``           | ``xx[a]``               |
| $x_{ab}^2 = \left(x_a-x_b\right)^2$   | ``xx(a, b)``        | ``xx[a,b]``             |
| $\frac{\partial f}{\partial x_a^i}$   | ``DD(f, (a, i))``   | ``DD[f,{a,i}]``         |

$x^i_a$ represents the set of $\mathbb{R}^d$ points on which the function depends. The first argument labels the point itself, while the second argument is the coordinate label. In python, both of them should be represented by SymPy symbols. Any vector parameters should also be encoded in this fashion.

The ``xx`` function is used to denote distances between points or between each point and the origin. Other notations include ``d`` for the dimension of the space and ``delta(i,j)``/``\[Delta][i,j]`` for the Kronecker delta.

Function ``DD`` computes the derivative of ``f`` with respect to $x^\mu_i$. For example,

```python
In: rdd.DD(xx(a, b) * x(b, i) + x(a, i), (a, j))
Out: 2 * (x(a, j) - x(b, j)) * x(b, i) + delta(j, i)
```
```mathematica
In: xx[a,b] x[b,i] + x[a,i] // DD[#,{a,j}] &
Out: 2 x[b,i] (x[a,j] - x[b,j]) + \[Delta][i,j]
```

For large-scale calculations in Mathematica, one might use the ``ParallelDD`` function instead.

This package also enables users to contract tensorial expressions over index pairs using the ``contract`` Python function or the ``Contract``/``ParallelContract`` mathematica functions, as demonstrated in the following example,

```python
In: rdd.contract(delta(i, k) * x(a, j) + x(a, i) * x(a, j) * x(b, k), (i, j))
Out: x(a, k) + x(b, k) * xx(a)
```
```mathematica
In: \[Delta][i,k] x[a,j] + x[a,i] x[a,j] x[b,k] // Contract[#,{i,j}] &
Out: x[a,k] + x[b,k] xx[a]
```
