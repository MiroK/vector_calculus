from parameter_domain import ParameterDomain
from vector_calculus.containers import Vector, Tensor
from vector_calculus.operators import dot, cross
from sympy import Number, symbols, diff, Matrix, sqrt
from numpy import array, ndarray
from numpy.linalg import det


# Symbols in terms of which the mapping is defined
__symbols__ = symbols('s, t, r')


class ParametrizedSet(object):
    '''Set described by mapping parameters from their ParameterDomain to R^d.'''

    def __init__(self, domain, mapping, orientation='+'):
        '''
        Construct the set. Orientation is used to construct normal if relevant.
        Positive for 2d-3d means normal = cross(d(mapping)/ds, d(mapping)/dt)
        Positve for 1d-2d means normal = rotate tangent counter-clockwie
        '''
        assert isinstance(domain, ParameterDomain)
        # Topological dimension
        self._tdim = len(domain) 

        # Map is a an iterable with length = gdim that specifies [x, y, z] as
        # functions of parameters of domain
        if not hasattr(mapping, '__len__'):
            mapping = (mapping, )

        # Check that the component uses known params or numbers
        for comp in mapping:
            comp_extras = comp.atoms() - domain.parameters
            
            params_only = len(comp_extras) == 0 
            extras_are_numbers = all(isinstance(item, (float, int, Number))
                                     for item in comp_extras)
            
            assert params_only or extras_are_numbers, 'Invalid mapping'

        # Geometrical dimension
        self._gdim = len(mapping)

        # I am only interested in at most 3d gdim
        assert 0 < self._gdim < 4, 'Geometrical dimension %s not supported' % self._gdim
        assert 0 < self._tdim < 4, 'Topological dimension %s not supported' % self._tdim
        assert self._tdim <= self._gdim, 'Topolog. dim > geometric. dim not allowed' 

        params = domain.parameters
        # Every mapping has a Jacobian but not every has normal and tangent
        self._n = None
        self._tau = None
       
        # Volumes, Square matrix only Jacobian 
        if self._tdim == self._gdim:
            Jac = Matrix([[diff(comp, var) for var in params] for comp in mapping])
            self._J = abs(Jac.det())
        # Curves and surfaces have normals or tangents in addition to Jacobian
        else:
            # Curves
            if self._tdim == 1:
                # Tagent
                self._tau = Vector([diff(comp, list(params)[0]) for comp in mapping])
                # Jacobian is length of tangent
                self._J = sqrt(sum(v**2 for v in self._tau))
                
                # And in 2d we can define a normal
                if self._gdim == 2:
                    R = Tensor([[0, -1], [1, 0]])
                    R = R if orientation == '+' else -R
                    self._n = dot(R, self._tau)

            # Surface in 3d has normal
            elif self._tdim == 2 and self._gdim == 3:
                u0 = Vector([diff(comp, list(params)[0]) for comp in mapping])
                u1 = Vector([diff(comp, list(params)[1]) for comp in mapping])

                self._n = cross(u0, u1)
                self._n = self._n if orientation == '+' else -self._n
                self._J = sqrt(sum(v**2 for v in self._n))

    @property
    def tdim(self):
        '''Topological dimension of set.'''
        return self._tdim

    @property
    def gdim(self):
        '''Geometrical dimension of set.'''
        return self._gdim

    @property
    def J(self):
        '''Jacobian.'''
        return self._J
    
    @property
    def tau(self):
        if self._tau is not None:
            return self._tau
        else:
            raise ValueError('No tangent for set with tdim=%d and gdim=%d' %\
                             (self.tdim, self.gdim))

    @property
    def n(self):
        if self._n is not None:
            return self._n
        else:
            raise ValueError('No normal for set with tdim=%d and gdim=%d' %\
                             (self.tdim, self.gdim))


class SimplexSet(ParametrizedSet):
    '''Simplex domain = smallest convex envelope of vertices.'''
    def __init__(self, vertices):
        assert 1 < len(vertices) < 5, 'Only line, triangle, tetrahedron'
        gdim = len(vertices[0])
        assert all(gdim == len(v) for v in vertices[1:]), \
            'Vertices of different length'

        # Check degeneracy
        mat = array([vertex if isinstance(vertex, ndarray) else array(vertex)
                    for vertex in vertices])
        mat -= mat[0, :]
        mat = mat[1:]
        assert det(mat.dot(mat.T)) > 1E-15, 'Degenerate simplex'

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


class CartesianSet(ParametrizedSet):
    '''Domain as cartesian product of intervals.'''
    def __init__(self, intervals):
        assert 0 < len(intervals) < 4, 'Only 1d, 2d, 3d'
        assert all(len(interval) == 2 for interval in intervals),\
            'Too many points in some interval'
        assert all(map(lambda pair: abs(pair[1] - pair[0]) > 1E-15, intervals)),\
            'Some interval is degenerate'
        assert all(map(lambda pair: pair[1] > pair[0], intervals)),\
            'Not increasing interval'

        tdim = len(intervals)
        smbls = __symbols__
        # Build domain
        domain = tuple((smbl, (-1, 1)) for smbl in smbls[:tdim])
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
    A = [0, 0, 0]
    B = [1, 0, 0]
    C = [0, 1, 0]
    D = [0, 0, 1]

    line = Line(A[:2], B[:2])
    tri = Triangle(A[:], B[:], C[:])
    tet = Tetrahedron(A, B, C, D)

    intrval = Interval(0, 1)
    # rect = Rectangle([0, 1], [0, 1])
    box = Box([0, 1], [0, 1], [0, 1])

    print box.J
    print line.J, line.tau, line.n
    print tri.J, tri.n

    from sympy import sin, cos, pi
    s = symbols('s')
    arc = ParametrizedSet(ParameterDomain((s, (0, pi))), (sin(s), cos(s)))
    print arc.J, arc.tau
