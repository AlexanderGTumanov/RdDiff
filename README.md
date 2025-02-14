# Coordinate differentiation for Mathematica
A Wolfram Mathematica implementation of coordinate differentiation in d-dimensional falt spaces. 
## Installation
Clone this repository using
```console
git clone https://github.com/AlexanderGTumanov/coordinate-differentiation-for-mathematica.git
```
This will create a new directory ``coordinate-differentiation-for-mathematica-for-mathematica`` in the current working directory. In your notebook file, add
```mathematica
SetDirectory["location_of_the_coordinate-differentiation-for-mathematica_folder"];
<<CoordinateDifferentiation`;
```
## Usage
This package is designed to perform various operations on functions, vectors, and tensors in a dd-dimensional flat space. The metric signature is irrelevant, as the code does not rely on it. That being said, it is the user's responsibility to ensure the correct placement of upper and lower indices. The package's notation can be summarized in the following table,

| Object                                | Notation            |
| --------                              | -------             |
| $x^\mu_i$                             | ``x[i,\[Mu]]``      |
| $x_i^2 = x_i^\mu\ x_{i\mu}$           | ``xx[i]``           |
| $x_{ij}^2 = \left(x_i-x_j\right)^2$   | ``xx[i,j]``         |
| $\frac{\partial f}{\partial x_i^\mu}$ | ``DD[f,{i,\[Mu]}]`` |

``x[i,\[Mu]]`` represents the set of variables on which the function depends. The first argument labels the variable itself—it does not have to be a number but must be unique for each variable. Any vector parameters should also be included in this set. The package does not support parameters with more complex tensor structures.

The ``xx`` function is used to denote distances between points or between each point and the origin. Other notations include ``d`` for the dimension of the space and ``\[Delta][\[Mu],\[Nu]]`` for the Kronecker delta.

Function ``DD[f,{i,\[Mu]}]`` computes the derivative of ``f`` with respect to $x^\mu_i$. For example,

```mathematica
In: xx[i,j] x[j,\[Mu]] + x[i,\[Mu]] // DD[#,{i,\[Nu]}] &
Out: 2 x[j,\[Mu]] (x[i,\[Nu]] - x[j,\[Nu]]) + \[Delta][\[Mu],\[Nu]]
```

For large-scale calculations, one might use the ``ParallelDD`` function instead.

This package also allows the user to contract tensorial expressions over pairs of indices. To do this, one must identify the pair of indices within the expression and apply the ``Contract`` (``ParallelContract``) function, as demonstrated in the following example,

```mathematica
In: \[Delta][\[Mu],\[Alpha]] x[i,\[Nu]] + x[i,\[Mu]] x[i,\[Nu]] x[j,\[Alpha]] /. \[Nu] -> \[Mu] // Contract[#,\[Mu]] &
Out: x[i,\[Alpha]] + x[j,\[Alpha]] xx[i]
```

This function can also contract multiple pairs of indices simultaneously if the second argument is provided as a list.
