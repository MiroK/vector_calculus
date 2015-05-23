from parameter_domain import ParameterDomain
from sympy import Number, symbols


# Symbols in terms of which the mapping is defined
__symbols__ = symbols('s, t, r')


class ParametrizedSet(object):
    '''Set described by mapping parameters from their ParameterDomain to R^d.'''

    def __init__(self, domain, mapping):
        '''Construct the set.'''
        assert isinstance(domain, ParameterDomain)
        # Topological dimension
        self._tdim = len(domain) 

        # Map is a an iterable with length = gdim that specifies [x, y, z] as
        # functions of parameters of domain
        if not hasattr(mapping, '__len__'):
            mapping = (mapping, )

        # Check that the component uses known params or numbers
        for comp in mapping:
            comp_extras = comp.atoms() - domain.parameters()
            
            params_only = len(comp_extras) == 0 
            extras_are_numbers = all(isinstance(item, (float, int, Number))
                                     for item in comp_extras)
            
            assert params_only or extras_are_numbers, 'Invalid mapping'

        # Geometrical dimension
        self._gdim = len(mapping)

    @property
    def tdim(self):
        '''Topological dimension of set.'''
        return self._tdim

    @property
    def gdim(self):
        '''Geometrical dimension of set.'''
        return self._gdim

    #FIXME, eval, call or something like that, something to be used with subs
    #FIXME, Jacobian, normal, ...


## Simplex
class SimplexSet(ParametrizedSet):
    '''Simplex domain = smallest convex envelope of vertices.'''
    def __init__(self, vertices):
        assert 1 < len(vertices) < 5, 'Only line, triangle, tetrahedron'
        gdim = len(vertices[0])
        assert all(gdim == len(v) for v in vertices[1:]), \
            'Vertices of different length'

        smbls = __symbols__
        tdim = len(vertices) - 1
        # Build mapping, A*(1-s) + B*s etc
        foo = [1-sum([s for s in smbls[:tdim]])] + list(smbls)
        mapping = tuple(sum(vtx[dim]*sym for vtx, sym in zip(vertices, foo))
                        for dim in range(gdim))

        # Build domain (s, (0, 1)), (r, (0, 1-s)), ...
        domain = tuple((smbls[dim], (0, 1-sum(s for s in smbls[:dim])))
                       for dim in range(tdim))
        domain = ParameterDomain(*domain)

        #FIXME, check degeneracy using Jacobian which is constant

        ParametrizedSet.__init__(self, domain, mapping)


class Line(SimplexSet):
    '''
    Line between points A, B in R^d characterized by x=A*s + (B-A)*(1-s) for s
    in [0, 1]. Here x is x, (x, y) or (x, y, z) for d=1, 2, 3. 
    '''
    def __init__(self, A, B):
        SimplexSet.__init__(self, [A, B])


class Triangle(SimplexSet):
    '''
    Triangle formed by vertices A, B, C in R^d, d=2, 3. The map is
    A(1-s-t) + Bs + Ct for s in [0, 1] and t in [0, 1-s].
    '''
    def __init__(self, A, B, C):
        SimplexSet.__init__(self, [A, B, C])


class Tetrahedron(SimplexSet):
    '''
    Tetrahedron formed by vertices A, B, C, D in R^3. The map is
    A(1-s-t-r) + Bs + Ct + Dr for s in [0, 1], t in [0, 1-s] and r in [0, 1-s-t].
    '''
    def __init__(self, A, B, C, D):
        SimplexSet.__init__(self, [A, B, C, D])


## Cartesian
class CartesianSet(ParametrizedSet):
    '''Domain as cartesian product of intervals.'''
    def __init__(self, intervals):
        assert 0 < len(intervals) < 4, 'Only 1d, 2d, 3d'
        assert all(len(interval) == 2 for interval in intervals),\
            'Too many points in some interval'
        assert all(map(lambda pair: abs(pair[1] - pair[0]) > 1E-13, intervals)),\
            'Some interval is degenerate'
        assert all(map(lambda pair: pair[1] > pair[0], intervals)),\
            'Not increasing interval'

        tdim = len(intervals)
        smbls = __symbols__
        # Build domain
        domain = tuple((smbl, (-1, 1)) for smbl in smbls)
        domain = ParameterDomain(*domain)

        gdim = tdim
        # Build mapping
        xi_map = lambda (interval, sym): 0.5*(interval[0]*(1-sym) + interval[1]*(1+sym))
        mapping = tuple(map(xi_map, zip(intervals, smbls)))

        ParametrizedSet.__init__(self, domain, mapping)


class Interval(CartesianSet):
    '''Interval [a, b] desribed as x = 0.5*a(1-s) + 0.5*b(1+s), s in [-1, 1].'''
    def __init__(self, a, b):
        CartesianSet.__init__(self, [[a, b]])


class Rectangle(CartesianSet):
    '''
    Rectanle [a0, b0] x [a1, b1] as x = 0.5*a0*(1-s) + 0.5*b0*(1+s) and
    y = 0.5*a1*(1-t) + 0.5*b1*(1+t) for (s, t) in [-1, 1] x [-1, 1].
    '''
    def __init__(self, xI, yI):
        CartesianSet.__init__(self, [xI, yI])


class Box(CartesianSet):
    '''
    Rectanle [a0, b0] x [a1, b1] x [a2, b2] as x = 0.5*a0*(1-s) + 0.5*b0*(1+s) and
    y = 0.5*a1*(1-t) + 0.5*b1*(1+t) and z = 0.5*a2*(1-t) + 0.5*b2*(1+t) for 
    (s, t, r) in [-1, 1] x [-1, 1] x [-1, 1].
    '''
    def __init__(self, xI, yI, zI):
        CartesianSet.__init__(self, [xI, yI, zI])


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    A = [0, 0, -1]
    B = [1, 1, 2]
    C = [1, 2, 3]
    D = [3, 5, 2]

    line = Line(A[:1], B[:1])
    tri = Triangle(A, B, C)
    tet = Tetrahedron(A, B, C, D)

    intrval = Interval(0, 1)
    rect = Rectangle([0, 1], [0, 1])
    box = Box([0, 1], [-1, 1], [0, 1])
