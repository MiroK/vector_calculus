from measure import Measure


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

    print 1*VolumeMeasure(Interval(1, 10))
    print quad(lambdify(x, S(1)), [1, 10])
