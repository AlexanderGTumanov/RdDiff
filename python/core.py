import sympy as sp

# d: dimension, x[i, a]: a-th component of vector x_i, xx[i, j]: (x_i - x_j)^2

d = sp.symbols('d')
x = sp.Function('x')
xx = sp.Function('xx', commutative = True)
delta = sp.Function('delta', commutative = True)

d.__doc__ = """ This symbol is reserved for the dimension of space."""
x.__doc__ = """ This function represents the a-th component of vector x_i."""
xx.__doc__ = """ This function represents distances between points. xx(a, b) is the distance between x_a and x_b, xx(a) is the distance between x_a and the origin (length of x_a)"""
delta.__doc__ = """ This function represents the Kronecker delta of two indexes."""

def arg_complement(obj, arg):
    """
    This function is used internally to extract the second argument of x, xx, and delta functions, when given the first.
    """
    arg_list = list(obj.atoms())
    return arg_list[1] if arg_list[0] == arg else arg_list[0]

def DD(expr, dir):
    """
    This function takes the derivative of an expression with respect to a coordiante vector.

    Parameters:
    expr (sp.Basic): The SymPy expression to be differentiated.
    dir (tuple of two sp.Symbol): A tuple/list with exactly two SymPy symbols: the first is the label (a) of x_a^i, the second is its component (i)

    Returns:
    sp.Basic
    """

    # listing all xx's with dir[0] as one of the arguments
    all_xx_instances = {atom for atom in expr.atoms(sp.Function) if atom.func == xx}
    xx_instances = {match for match in all_xx_instances if len(match.args) == 2 and (match.args[0] == dir[0] or match.args[1] == dir[0])}

    # listing all x's with dir[0] as the first argument
    all_x_instances = {atom for atom in expr.atoms(sp.Function) if atom.func == x}
    x_instances = {match for match in all_x_instances if len(match.args) == 2 and match.args[0] == dir[0]}
    
    # adding up all contributions from differenting xx's
    xx_part = 2 * sp.diff(expr, xx(dir[0])) * x(dir[0], dir[1])
    for el in xx_instances:
        xx_part += 2 * sp.diff(expr, el) * (x(dir[0], dir[1]) - x(arg_complement(el, dir[0]), dir[1]))
    
    # adding up all contributions from differenting x's
    x_part = 0
    for el in x_instances:
        x_part += expr.expand().coeff(el, 1) * delta(dir[1], arg_complement(el, dir[0]))
    
    return xx_part + x_part

def contract(expr, ind):
    """
    This function takes the derivative of an expression with respect to a coordiante vector.

    Parameters:
    expr (sp.Basic): The SymPy expression to be differentiated.
    dir (tuple of two sp.Symbol): A tuple/list with exactly two SymPy symbols: the first is the label (a) of x_a^i, the second is its component (i)

    Returns:
    sp.Basic
    """

    # listing all x's with ind[0] and ind[1] as first arguments
    all_x_instances = {atom for atom in expr.atoms(sp.Function) if atom.func == x}
    x0_instances = {match for match in all_x_instances if len(match.args) == 2 and match.args[1] == ind[0]}
    x1_instances = {match for match in all_x_instances if len(match.args) == 2 and match.args[1] == ind[1]}

    # listing all delta's with ind[0] as one of the arguments (except delta(ind[0], ind[1])
    all_delta_instances = {atom for atom in expr.atoms(sp.Function) if atom.func == delta}
    all_delta_instances.discard(delta(ind[0], ind[1]))
    all_delta_instances.discard(delta(ind[1], ind[0]))
    delta0_instances = {match for match in all_delta_instances if len(match.args) == 2 and (match.args[0] == ind[0] or match.args[0] == ind[0])}
    delta1_instances = {match for match in all_delta_instances if len(match.args) == 2 and (match.args[0] == ind[1] or match.args[0] == ind[1])}

    # adding up all contributions from products of x_a^ind[0] x_b^ind[1]-type terms
    x_x_part = 0
    for el0 in x0_instances:
        for el1 in x1_instances:
            if arg_complement(el0, ind[0]) == arg_complement(el1, ind[1]):
                x_x_part += xx(arg_complement(el0, ind[0])) * expr.expand().coeff(el0, 1).expand().coeff(el1, 1)
            else:
                x_x_part += sp.Rational(1, 2) * (xx(arg_complement(el0, ind[0])) + xx(arg_complement(el1, ind[1])) - xx(arg_complement(el0, ind[0]), arg_complement(el1, ind[1]))) * expr.expand().coeff(el0, 1).expand().coeff(el1, 1)

    # adding up all contributions from products of x_a^ind[0] delta_i^ind[1]-type terms
    x_delta_part = 0
    for el0 in x0_instances:
        for el1 in delta1_instances:
            x_delta_part += x(arg_complement(el0, ind[0]), arg_complement(el1, ind[1])) * expr.expand().coeff(el0, 1).expand().coeff(el1, 1)

    # adding up all contributions from products of delta_i^ind[0] x_a^ind[1]-type terms
    delta_x_part = 0
    for el0 in delta0_instances:
        for el1 in x1_instances:
            delta_x_part += x(arg_complement(el1, ind[1]), arg_complement(el0, ind[0])) * expr.expand().coeff(el0, 1).expand().coeff(el1, 1)

    # adding up all contributions from products of delta_i^ind[0] delta_j^ind[1]-type terms
    delta_delta_part = 0
    for el0 in delta0_instances:
        for el1 in delta1_instances:
            delta_delta_part += delta(arg_complement(el0, ind[0]), arg_complement(el1, ind[1])) * expr.expand().coeff(el0, 1).expand().coeff(el1, 1)
    
    return x_x_part + x_delta_part + delta_x_part + delta_delta_part + d * expr.expand().coeff(delta(ind[0], ind[1])) + d * expr.expand().coeff(delta(ind[1], ind[0]))