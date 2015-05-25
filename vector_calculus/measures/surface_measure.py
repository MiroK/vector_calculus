from sympy import Expr, Number, NumberSymbol
from vector_calculus.containers import Vector, Tensor
from vector_calculus.operators import dot, inner
from measure import Measure


class SurfaceMeasure(Measure):
    '''
    SurfaceMeasure is defined over sets which have same geometrical dimension
    one higher than topological dimension.
    '''
    
    def __init__(self, domain):
        assert domain.gdim == domain.tdim + 1, \
            'Invalid domain tdim(%d) + 1!= gdim(%d)' % (domain.tdim, domain.gdim)
        Measure.__init__(self, domain)


    def __rmul__(self, integrand):
        '''Integrate over domain.'''
        # Scalar integral is f(x(s, t), y(s, t))*|d(x, y)/ds x d(x, y)/dt| ds dt
        # Result is number
        if isinstance(integrand, (Expr, Number, NumberSymbol, int, float)):
            # Note that Jacobian for thich domain is defines as size of the
            # normal vector
            integrand = integrand*self.domain.J
            return self(integrand)
        # Vector integral is 
        # inner(f(x(s, t), y(s, t)), (d(x, y)/ds x d(x, y)/dt)) ds dt
        # Result is number
        elif isinstance(integrand, Vector):
            assert len(integrand) == self.domain.gdim, 'Gdim mismatch'
            integrand = inner(integrand, self.domain.n)
            return self(integrand)
        # Tensor integral is
        # dot(f(x(s, t), y(s, t)), (d(x, y)/ds x d(x, y)/dt)) ds dt
        # Result is vector
        elif isinstance(integrand, Tensor):
            assert len(integrand) == self.domain.gdim, 'Gdim mismatch'
            integrand = dot(integrand, self.domain.n)
            return Vector([vi*self for vi in integrand])
        # Nope 
        else:
            raise TypeError('No surface integral of type %s' % type(integrand))


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    from parametrized_set import Triangle
    from sympy import symbols

    A = [1, 0, 0]
    B = [0, 1, 0]
    C = [0, 0, 1]
    dS = SurfaceMeasure(Triangle(A, B, C))

    x, y, z = symbols('x, y, z')
    f = x + y + z
    v = Vector([x, y, z])

    print v*dS

    T = Tensor([[x, 0, 0], [0, y, 0], [0, 0, z]])
    print T*dS
