from copy import copy
from sympy import Matrix


class Vector(object):
    '''Symbolic vector.'''

    def __init__(self, block):
        '''Vector is a fancy list'''
        assert len(block) == 2 or len(block) == 3, 'Only 2d and 3d vectors'
        self.u = copy(block)

    def __getitem__(self, i):
        '''Extract component.'''
        return self.u[i]

    def __len__(self):
        '''Length of vector.'''
        return len(self.u)

    def __str__(self):
        '''String representation.'''
        return self.u.__str__()

    def __add__(self, v):
        '''Add two vectors.'''
        return Vector([ui+vi for ui, vi in zip(self, v)])

    def __sub__(self, v):
        '''Subtract two vectors.'''
        return Vector([ui-vi for ui, vi in zip(self, v)])

    def __mul__(self, a):
        '''Multiply by scalar.'''
        return Vector([ui*a for ui in self])

    def __rmul__(self, a):
        '''Multiply by scalar.'''
        return Vector([ui*a for ui in self])

    def __div__(self, a):
        '''Divide by scalar.'''
        return Vector([ui/a for ui in self])

    def __neg__(self):
        '''Multiply by -1.'''
        return self*-1

    def __eq__(self, v):
        '''Check equality. Depends on == in sympy so use with caution.'''
        return all(ui == vi for ui, vi in zip(self, v))

    def subs(self, values):
        '''Substitute each component.'''
        return Vector([ui.subs(values) for ui in self])
    
    # FIXME: Maybe start with Matrix
    def as_matrix(self):
        '''Return copy as sympy Matrix.'''
        return Matrix(self.u)


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    from sympy import symbols

    print 'TODO: TESTS!'
    x = symbols('x')

    u = Vector([x, x+1])
    print u.u
