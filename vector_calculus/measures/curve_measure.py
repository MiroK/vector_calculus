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

    def __init__(self, *domain):
        if len(domain) == 2:
            # dL(A, B);  line between A, B
            if hasattr(domain[-1], '__len__'):
                A, B = domain[0], domain[1]
                CurveMeasure.__init__(self, Line(A, B))
            else:
                assert isinstance(domain[-1], int)
                # dL([[[], []], index); edge of rectangle
                if len(domain[0]) == 2:
                    [[a0, b0], [a1, b1]] = domain[0]
                    A = [a0, a1]
                    B = [b0, a1]
                    C = [b0, b1]
                    D = [a0, b1]

                    index = domain[-1]
                    assert 0 <= index < 4, 'Invalid index %d for rectangle' % index
                    # Edges are 0: A -> B, 1: B -> C, 2: C -> D, 3: D -> A 
                    CurveMeasure.__init__(self, Line(*{0: (A, B),
                                                       1: (B, C),
                                                       2: (C, D),
                                                       3: (D, A)}[index]))
                # dL([[], [], []], index); edge of box
                elif len(domain[0]) == 3:
                    [[a0, b0], [a1, b1], [a2, b2]] = domain[0]
                    A = [a0, a1, a2]
                    B = [b0, a1, a2]
                    C = [b0, b1, a2]
                    D = [a0, b1, a2]
                    E = [a0, a1, b2]
                    F = [b0, a1, b2]
                    G = [b0, b1, b2]
                    H = [a0, b1, b2]

                    index = domain[-1]
                    assert 0 <= index < 12, 'Invalid index %d for box' % index
                    CurveMeasure.__init__(self, Line(*{0: (A, B), 1: (B, C),
                                                       2: (C, D), 3: (D, A),
                                                       4: (A, E), 5: (B, F),
                                                       6: (C, G), 7: (D, H),
                                                       8: (E, F), 9: (F, G),
                                                       10: (G, H), 11: (H, E)}[index]))
                else:
                    raise ValueError('Invalid domain')
        # dL(A, B, C, index); edge of triangle
        elif len(domain) == 4 and isinstance(domain[-1], int):
            index = domain[-1]
            assert 0 <= index < 3, 'Invalid index %d for triangle' % index
            # Edges are 0: A -> B, 1: B -> C, 2: C -> A
            A, B, C = domain[0], domain[1], domain[2]
            CurveMeasure.__init__(self, Line(*{0: (A, B),
                                               1: (B, C),
                                               2: (C, A)}[index]))

        # dL(A, B, C, D, index); edge of tet
        elif len(domain) == 5 and isinstance(domain[-1], int):
            index = domain[-1]
            assert 0 <= index < 6, 'Invalid index %d for tetrahedron' % index
            A, B, C, D = domain[0], domain[1], domain[2], domain[3]
            # Edges are 0: A -> B, 1: B -> C, 2: C -> D, 3: D -> A, 4: B -> D,
            # 5: C -> A
            CurveMeasure.__init__(self, Line(*{0: (A, B),
                                               1: (B, C),
                                               2: (C, D),
                                               3: (D, A),
                                               4: (B, D),
                                               5: (C, A)}[index]))
        else:
            raise ValueError('Invalid domain')


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    from parametrized_set import Line
    from sympy import symbols

    A = [0, 0]
    B = [1, 0]
    C = [1, 1]


    ds = sum([dL(A, B), dL(B, C)], dL(C, A))

    print '>>', 1*ds
