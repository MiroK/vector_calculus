from vector import Vector
from sympy import Number
from copy import copy
from sympy import Matrix
from sympy import Number, NumberSymbol, Expr


class Tensor(object):
    '''Symbolic rank-2 tensor.'''

    def __init__(self, blocks):
        '''Build tensor from list or list or Vectors'''
        assert len(blocks) == 2 or len(blocks) == 3,\
            'Only 2d and 3d tensor allowed'
        dim = len(blocks)
        assert all(len(block) == dim for block in blocks),\
            'Row length does not match dim = %d' % dim
        # Finally build
        self.A = [block if isinstance(block, Vector) else Vector(block)
                  for block in blocks]

    def __getitem__(self, i):
        '''Extract component.'''
        return self.A[i]

    def __len__(self):
        '''Length of vector.'''
        return len(self.A)

    def __str__(self):
        '''String representation.'''
        return '[' + ' '.join([u.__str__() for u in self]) + ']'

    def __add__(self, B):
        '''Add two vectors.'''
        return Tensor([Ai+Bi for Ai, Bi in zip(self, B)])

    def __sub__(self, B):
        '''Subtract two vectors.'''
        return Tensor([Ai-Bi for Ai, Bi in zip(self, B)])

    def __mul__(self, a):
        '''Multiply by scalar.'''
        # Scalar
        if isinstance(a, (float, int, Number, NumberSymbol, Expr)):
            return Tensor([Ai*a for Ai in self])
        # Multiply two tensor
        elif isinstance(a, Tensor):
            assert len(a) == len(self),\
                'Incompatible tensor lenghts %d and %d' % (len(self), len(a))
            n = len(a)
            blocks = []
            for i in range(n):
                row = []
                for j in range(n):
                    row.append(sum([self[i][k]*a[k][j] for k in range(n)]))
                blocks.append(row)
            return Tensor(blocks)
        # No other
        else:
            return NotImplemented

    def __pow__(self, n):
        '''Operation for A**n.'''
        assert isinstance(n, int)
        assert n >= 0
        if n == 0:
            return copy(self)
        # Forget memory efficiency
        else:
            return reduce(lambda A, B: A*B, [self]*n)
        
    def __rmul__(self, a):
        '''Multiply by scalar.'''
        return Tensor([Ai*a for Ai in self])

    def __div__(self, a):
        '''Divide by scalar.'''
        return Tensor([Ai*(1./a) for Ai in self])

    def __neg__(self):
        '''Multiply by -1.'''
        return self*-1

    def __eq__(self, B):
        '''Check equality. Depends on == in sympy so use with caution.'''
        return all(Ai == Bi for Ai, Bi in zip(self, B))

    def subs(self, values):
        '''Substitute each component.'''
        return Tensor([Ai.subs(values) for Ai in self])

    # FIXME: Matrix is probably a better container
    def as_matrix(self):
        '''Return copy as sympy Matrix.'''
        return Matrix(self.A)
