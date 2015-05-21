from sympy import integrate, symbols
from copy import copy

xyz = symbols('x, y, z')


class dx(object):
    '''Volume integral in 2d, 3d.'''

    def __init__(self, *args):
        '''
        Constructors:

        dx([a0, b0], [a1, b1]) and dx([a0, b0], [a1, b1], [a2, b2]) are volume
        integrals over rectangle and box.

        dx([A, B, C]) and dx([A, B, C, D]) are volume integrals over triangle
        and tetrahedron.
        '''
        # Rectangular domains
        if len(args) == 2 or len(args) == 3:
            assert all(len(arg) == 2 for arg in args)
            dim = len(args)

            # 2d integration assumes integrand of x, y
            if dim == 2:
                self._dim = dim
                x, y = xyz[:dim]
                dx, dy = args
                self._integrate = \
                    lambda f: integrate(integrate(f, (x, dx[0], dx[1])),
                                        (y, dy[0], dy[1]))
            # 3d integration 
            else:
                self._dim = dim
                x, y, z = xyz[:dim]
                dx, dy, dz = args
                self._integrate = \
                    lambda f: integrate(integrate(integrate(f, (x, dx[0], dx[1])),
                                                  (y, dy[0], dy[1])),
                                        (z, dz[0], dz[1]))

    def __rmul__(self, integrand):
        '''Integrate integrand with the measure'''
        return self._integrate(integrand)
    
    @property
    def dim(self):
        '''Geometric dim of domain where the volume ingral is taken.'''
        return self._dim


    def __add__(self, other):
        '''Add two volume measures.'''
        assert isinstance(other, dx)
        assert self.dim == other.dim
        dxx = copy(self)
        dxx._integrate = lambda f: self._integrate(f) + other._integrate(f)
        return dxx


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    x, y = xyz[:2]
    f = (x**2+y**2)
    
    dX = dx([-1, 1], [-1, 1])
    dXX = dx([-1, 1], [-1, 2])

    print float(f*(dX+dXX))

    from sympy.mpmath import quad
    from sympy import lambdify

    print quad(lambdify([x, y], f), [-1, 1], [-1, 1]) + quad(lambdify([x, y], f), [-1, 1], [-1, 2])
