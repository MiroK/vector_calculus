from sympy import Expr, Number, NumberSymbol
from vector_calculus.containers import Vector, Tensor
from vector_calculus.operators import dot, inner
from measure import Measure


class CurveMeasure(Measure):
    '''
    CurveMeasure is defined over sets whose topological dimension equals 1.
    '''
    
    def __init__(self, domain):
        assert domain.tdim == 1, \
            'Invalid domain tdim(%d) != 1' % domain.tdim
        Measure.__init__(self, domain)


    def __rmul__(self, integrand):
        '''Integrate over domain.'''
        # Scalar integral is f(x(s), y(s))*|d(x, y)/ds| ds. Result is number
        if isinstance(integrand, (Expr, Number, NumberSymbol, int, float)):
            # Note that Jacobian for this domain is defines as size of the
            # tangent vector
            integrand = integrand*self.domain.J
            return self(integrand)
        # Vector ingral is inner(f(x(s), y(s)), d(x, y)/ds) ds. Result is number
        elif isinstance(integrand, Vector):
            assert len(integrand) == self.domain.gdim, 'Gdim mismatch'
            integrand = inner(integrand, self.domain.tau)
            return self(integrand)
        # Tensor ingral is dot(f(x(s), y(s)), d(x, y)/ds) ds. Result is vector
        elif isinstance(integrand, Tensor):
            assert len(integrand) == self.domain.gdim, 'Gdim mismatch'
            integrand = dot(integrand, self.domain.tau)
            return Vector([vi*self for vi in integrand])
        else:
            raise TypeError('No surface integral of type %s' % type(integrand))


class dL(CurveMeasure):
    '''
    Convenience function for defining line integrals over 'common' domains.
    '''
    # dL(A, B);  line between A, B
    # dL(A, B, C, index); edge of triangle
    # dL(A, B, C, D, index); edge of tet
    # dL([[[], []], index); edge of rectangle
    # dL([[], [], []], index); edge of box
    pass


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    from parametrized_set import Line
    from sympy import symbols

    A = [0, 0, 0]
    B = [1, 1, 1]

    x, y, z = symbols('x, y, z')
    f = Vector([x, y, z])

    dl = CurveMeasure(Line(A, B))
    print f*dl
