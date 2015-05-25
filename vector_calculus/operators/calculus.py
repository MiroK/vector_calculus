from vector_calculus.containers import *
from sympy import symbols, Expr, S
from linalg import tr, dot

# These are cannonical variables of cartesian coordinate system
xyz = symbols('x, y, z')


def div(u):
    '''Divergence of vector --> scalar. Divergence of tensor --> vector.'''
    # Vector
    if isinstance(u, Vector):
        return sum((Dx(ui, vari) for ui, vari in zip(u, xyz)), S(0))
    # Tensor, recurse rows
    elif isinstance(u, Tensor):
        return Vector([div(ui) for ui in u])
    else:
        raise TypeError('Only divergence of vector or tensor allowed.')


def grad(u, dim=None):
    '''Gradient of vector --> tensor. Gradient of scalar --> vector.'''
    # Scalar
    if isinstance(u, Expr):
        # Infer dim from arguments of u
        if dim is None:
            # If there is no z dependence this is most likely a 2d scalar
            if xyz[2] not in u.atoms():
                dim = 2
            else:
                dim = 3

        return Vector([Dx(u, var) for var in xyz[:dim]])
    # Vector, recurse rows
    elif isinstance(u, Vector):
        dim = len(u)
        return Tensor([grad(ui, dim) for ui in u])
    else:
        raise ValueError('Only gradient of scalar or vector allowed.')


def curl(u):
    '''Curl of 3d vector --> vector. Curl 2d vecror --> scalar.'''
    assert isinstance(u, Vector), 'Need vector for curl'
    if len(u) == 3:
        return -Vector([Dx(u[1], xyz[2]) - Dx(u[2], xyz[1]),
                        Dx(u[2], xyz[0]) - Dx(u[0], xyz[2]),
                        Dx(u[0], xyz[1]) - Dx(u[1], xyz[0])])
    else:
        R = Tensor([[0, 1], [-1, 0]])
        return div(dot(R, u))


def rot(u, orientation='+'):
    '''Rotation of 2d scalar --> 2d vector. Default is counter-clockwise rot.'''
    assert isinstance(u, Expr), 'Can only take rot of scalar'
    assert xyz[2] not in u.atoms(), 'Scalar must be function of x, y only'

    R = Tensor([[0, -1], [1, 0]])
    R = R if orientation == '+' else -R
    return dot(R, grad(u))


def Dx(u, var):
    '''Partial derivative of u w.r.t to var. Shape is maintained.'''
    if isinstance(u, (int, float)):
        u = S(u)
    # Scalar
    if isinstance(u, Expr):
        return u.diff(var, 1)
    # Vector
    elif isinstance(u, Vector):
        return Vector([Dx(ui, var) for ui in u])
    # Tensor
    elif isinstance(u, Tensor):
        return Tensor([Dx(ui, var) for ui in u])
    else:
        raise TypeError('Cannot take derivatieve of type %s' % type(u))
