from vector_calculus.containers import *
from sympy import Rational
from numpy import eye


def tr(A):
    'Trace of tensor.'
    return sum(A[i][i] for i in range(len(A)))


def transpose(A):
    'Return transpose of A.'
    n = len(A)
    blocks = [[] for i in range(n)]
    for row in A:
        for j in range(n):
            blocks[j].append(row[j])
    return Tensor(blocks)


def sym(A):
    'Return symmetrized tensor from A.'
    return Rational(1, 2)*(A+transpose(A))


def skew(A):
    'Return skew symmetrized tensor from A.'
    return Rational(1, 2)*(A-transpose(A))


def Id(n):
    'Idenity in R^n'
    return Tensor(eye(n))


def deviatoric(A):
    'Return deviatoric part of A.'
    n = len(A)
    return A - tr(A)*Id(n)/n


def commutator(A, B):
    'Commutator [A, B].'
    assert isinstance(A, Tensor) and isinstance(B, Tensor), 'Need two tensors'
    return A*B - B*A


def det(A):
    'Determinant of A.'
    return A.as_matrix().det()


def inv(A):
    'Formal inverse of A.'
    return A.as_matrix().inv()


def cross(u, v):
    'Cross product of two vectors --> vector.'
    assert isinstance(u, Vector) and isinstance(v, Vector), 'Need two vectors'
    assert len(u) == len(v), 'Need two vectors of same length'
    assert len(u) == 3, 'Need two vectors of lenght 3'

    return Vector([u[1]*v[2] - u[2]*v[1],
                   u[2]*v[0] - u[0]*v[2],
                   u[0]*v[1] - u[1]*v[2]])


def outer(u, v):
    'Outer product of two vectors --> tensor.'
    assert isinstance(u, Vector) and isinstance(v, Vector), 'Need two vectors'
    assert len(u) == len(v), 'Need two vectors of same length'
    
    n = len(u)
    blocks = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(u[i]*u[j])
        blocks.append(row)

    return Tensor(blocks)


def inner(u, v):
    'Inner product of two vectors or two tensors --> number.'
    assert all(isinstance(arg, Vector) for arg in (u, v)) or \
        all(isinstance(arg, Tensor) for arg in (u, v)),\
        'Arguments must be two vectors or two tensors'

    if isinstance(u, Vector):
        return sum((ui*vi for ui, vi in zip(u, v)))
    else:
        return tr(transpose(u)*v)


def dot(A, u):
    'Dot product between vector/tensor and vector/tensor.'
    assert isinstance(A, (Vector, Tensor)) and isinstance(A, (Vector, Tensor)),\
        'Dot product is between vector and tensors'
    
    # A.u
    if isinstance(A, Tensor):
        assert isinstance(u, Vector), 'A is Tensor but u in not a Vector'
        assert len(A) == len(u), 'Incompatible dimension'
        return Vector([inner(Ai, u) for Ai in A])
    # u.A is defined via transpose. We only have a row vector
    else:
        assert isinstance(u, Tensor), 'A is Vector but u in not a Tensor'
        return dot(transpose(u), A)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    from sympy import Symbol, simplify

    x = Symbol('x')
    A = deviatoric(Tensor([[10*x, 2, 3], [-1, -2, -3], [0, 0, 1]]))
    print tr(A)

    print A
    print transpose(A)
    print sym(A)
    print skew(A)
    print deviatoric(A)
    print commutator(A, A)
    print det(A)
    print inv(A)

    print simplify(inner(A, Id(len(A))))


    B = sym(A)
    v = Vector([1, 1, 1])
    print dot(B, v)
    print dot(v, B)
