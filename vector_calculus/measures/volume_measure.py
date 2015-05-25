from measure import Measure
from parametrized_set import Triangle, Tetrahedron, Rectangle, Box, Interval


class VolumeMeasure(Measure):
    '''
    VolumeMeasure is defined over sets which have same geometrical and
    topological dimension.
    '''
    
    def __init__(self, domain):
        assert domain.gdim == domain.tdim, \
            'Invalid domain tdim(%d) != gdim(%d)' % (domain.tdim, domain.gdim)
        Measure.__init__(self, domain)


    def __rmul__(self, integrand):
        '''Integrate over domain.'''
        # Add Jacobian
        integrand = integrand*self.domain.J
        return self(integrand)

        # FIXME convenience functions for Vector and Tensor


class dV(VolumeMeasure):
    '''
    Convenience function for defining volume integrals over 'common' domains.
    '''
    def __init__(self, *domain):
        # dV(A, B, C) over triangle
        if len(domain) == 3:
            VolumeMeasure.__init__(self, Triangle(*domain))
        # dV(A, B, C, D) over tetrahedron
        elif len(domain) == 4:
            VolumeMeasure.__init__(self, Tetrahedron(*domain))
        # Cartesian domain
        elif len(domain) == 1:
            domain = domain[0]
            # dV([[a0, b0]]) over interval
            if len(domain) == 1:
                VolumeMeasure.__init__(self, Interval(domain[0][0], domain[0][1]))
            # dV([[a0, b0], [a1, b1]]) over rectangle
            elif len(domain) == 2:
                VolumeMeasure.__init__(self, Rectangle(*domain))
            # dV([[a0, b0], [a1, b1], [a2, b2]]) over box
            elif len(domain) == 3:
                VolumeMeasure.__init__(self, Box(*domain))


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    from parametrized_set import Rectangle, Tetrahedron, Interval
    from sympy import lambdify, symbols, S
    from sympy.mpmath import quad

    x, y, z = symbols('x, y, z')
    f = x + y + z
    # print quad(lambdify([x, y], f), [0, 1], [1, 2])

    dx = VolumeMeasure(Rectangle([0, 1], [1, 3]))
    # print f*dx

    A = [0, 0, 0]
    B = [1, 0, 0]
    C = [0, 1, 0]
    D = [0, 0, 1]
    dx = VolumeMeasure(Tetrahedron(A, B, C, D))
    # print 1*dx
    # print 1*dV(A, B, C, D)

    f = x + y + z
    print f*dV([[0, 3], [0, 2], [0, 1]])
    f = x + y
    print f*dV([[0, 3], [0, 2]])
    f = x
    print f*dV([[1, 10]])

    f = x + y + z
    print quad(lambdify([x, y, z], f), [0, 3], [0, 2], [0, 1])
    f = x + y
    print quad(lambdify([x, y], f), [0, 3], [0, 2])
    f = x
    print quad(lambdify([x], f), [1, 10])

    # print 1*VolumeMeasure(Interval(1, 10))
    # print quad(lambdify(x, S(1)), [1, 10])
