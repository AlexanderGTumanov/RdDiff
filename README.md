# Coordinate differentiation for Mathematica
This package/library provides tools for performing differentiation and index contractions on tensor structures in $\mathbb{R}^d$. The structures are assumed to respect $O(d)$ symmetry and therefore depend only on coordinate vectors of different points and the distances between them. The Python library is built on top of SymPy.
## Installation

### Python
Clone this repository using
```console
git clone https://github.com/AlexanderGTumanov/coordinate-differentiation-for-mathematica.git
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

| Object                                | Python Notations    |  Mathematica Notations  |
| --------                              | -------             | -------                 |
| $x^i_a$                               | ``x(a,i)``          | ``x[a,i]``              |
| $x_a^2 = x_a^i\ x_{ai}$               | ``xx(a)``           | ``xx[a]``               |
| $x_{ab}^2 = \left(x_a-x_b\right)^2$   | ``xx(a,b)``         | ``xx[a,b]``             |
| $\frac{\partial f}{\partial x_a^i}$   | ``DD(f,(a,i))``     | ``DD[f,{a,i}]``         |

$x^i_a$ represents the set of $\mathbb{R}^d$ points on which the function depends. The first argument labels the point itself, while the second argument is the coordinate label. In python, both of them should be represented by SymPy symbols. Any vector parameters should also be encoded in this fashion.

The ``xx`` function is used to denote distances between points or between each point and the origin. Other notations include ``d`` for the dimension of the space and ``delta(i,j)``/``\[Delta][i,j]`` for the Kronecker delta.

Function ``DD`` computes the derivative of ``f`` with respect to $x^\mu_i$. For example,

```mathematica
In: xx[i,j] x[j,\[Mu]] + x[i,\[Mu]] // DD[#,{i,\[Nu]}] &
Out: 2 x[j,\[Mu]] (x[i,\[Nu]] - x[j,\[Nu]]) + \[Delta][\[Mu],\[Nu]]
```
```python
In: xx[i,j] x[j,\[Mu]] + x[i,\[Mu]] // DD[#,{i,\[Nu]}] &
Out: 2 x[j,\[Mu]] (x[i,\[Nu]] - x[j,\[Nu]]) + \[Delta][\[Mu],\[Nu]]
```

For large-scale calculations in Mathematica, one might use the ``ParallelDD`` function instead.

This package also allows the user to contract tensorial expressions over pairs of indices. To do this, one must identify the pair of indices within the expression and apply the ``Contract`` (``ParallelContract``) function, as demonstrated in the following example,

```mathematica
In: \[Delta][\[Mu],\[Alpha]] x[i,\[Nu]] + x[i,\[Mu]] x[i,\[Nu]] x[j,\[Alpha]] /. \[Nu] -> \[Mu] // Contract[#,\[Mu]] &
Out: x[i,\[Alpha]] + x[j,\[Alpha]] xx[i]
```

This function can also contract multiple pairs of indices simultaneously if the second argument is provided as a list.
