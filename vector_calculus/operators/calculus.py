from vector_calculus.containers import *
from sympy import symbols, Expr, S
from linalg import tr, dot

# These are cannonical variables of cartesian coordinate system
xyz = symbols('x, y, z')


def div(u):
    '''Divergence of vector --> scalar. Divergence of tensor --> vector.'''
    # Vector
    if isinstance(u, Vector):
        return sum((dx(ui, vari) for ui, vari in zip(u, xyz)), S(0))
    # Tensor, recurse rows
    elif isinstance(u, Tensor):
        return Vector([div(ui) for ui in u])
    else:
        raise TypeError('Only divergence of vector or tensor allowed.')


def grad(u):
    '''Gradient of vector --> tensor. Gradient of scalar --> vector.'''
    # Scalar
    if isinstance(u, Expr):
        # If there is no z dependence this is most likely a 2d scalar
        if xyz[2] not in u.atoms():
            dim = 2
        else:
            dim = 3

        return Vector([dx(u, var) for var in xyz[:dim]])
    # Vector, recurse rows
    elif isinstance(u, Vector):
        dim = len(u)
        return Tensor([grad(ui) for ui in u])
    else:
        raise ValueError('Only gradient of scalar or vector allowed.')


def curl(u):
    '''Curl of 3d vector --> vector. Curl 2d vecror --> scalar.'''
    assert isinstance(u, Vector), 'Need vector for curl'
    if len(u) == 3:
        return Vector([dx(u[1], xyz[2]) - dx(u[2], xyz[1]),
                       dx(u[2], xyz[0]) - dx(u[0], xyz[2]),
                       dx(u[0], xyz[1]) - dx(u[1], xyz[2])])
    else:
        R = Tensor([[0, 1], [-1, 0]])
        return div(dot(R, u))


def rot(u):
    '''Rotation of 2d scalar --> 2d vector.'''
    assert isinstance(u, Expr), 'Can only take rot of scalar'
    assert xyz[2] not in u.atoms(), 'Scalar must be function of x, y only'

    R = Tensor([[0, 1], [-1, 0]])
    return dot(R, grad(u))


def dx(u, var):
    '''Partial derivative of u w.r.t to var. Shape is maintained.'''
    # Scalar
    if isinstance(u, Expr):
        return u.diff(var, 1)
    # Vector
    elif isinstance(u, Vector):
        return Vector([dx(ui, var) for ui in u])
    # Tensor
    elif isinstance(u, Tensor):
        return Tensor([dx(ui, var) for ui in u])
    else:
        raise TypeError('Cannot take derivatieve of type %s' % type(u))


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    from sympy import symbols, sin, exp

    x, y, z = xyz

    u = Vector([x**2, y**3])
    print dx(u, x)

    f = 2*x**2 + y**3
    print f
    print grad(f)
    print grad(grad(f))

    print div(u)
    print div(grad(f))

    print curl(Vector([x, y, z]))
    print curl(Vector([y, x**2]))
    print rot(f)
